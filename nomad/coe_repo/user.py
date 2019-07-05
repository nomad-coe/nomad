# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Dict
from passlib.hash import bcrypt
from sqlalchemy import Column, Integer, String, ForeignKey, DateTime
from sqlalchemy.orm import relationship
import datetime
import jwt
import random
import string

from nomad import infrastructure, config, utils

from .base import Base


class Session(Base):  # type: ignore
    __tablename__ = 'sessions'

    token = Column(String)
    user_id = Column(String, ForeignKey('users.user_id'), primary_key=True)
    user = relationship('User')


class LoginException(Exception):
    """ Exception that is raised if the user could not be logged in despite present
    credentials. """
    pass


class Affiliation(Base):
    __tablename__ = 'affiliations'
    a_id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    address = Column(String)
    email_domain = Column(String)


class User(Base):  # type: ignore
    """
    SQLAlchemy model class that represents NOMAD-coe repository postgresdb *users*.
    Provides functions for authenticating via password or session token.

    It is not intended to create or update users. This should be done via the
    NOMAD-coe repository GUI.
    """
    __tablename__ = 'users'

    user_id = Column(Integer, primary_key=True)
    affiliation_id = Column(Integer, ForeignKey('affiliations.a_id'), name='affiliation')
    email = Column(String)
    first_name = Column(String, name='firstname')
    last_name = Column(String, name='lastname')
    affiliation = relationship('Affiliation', lazy='joined')
    password = Column(String)
    created = Column(DateTime)

    _token_chars = string.ascii_uppercase + string.ascii_lowercase + string.digits

    def __repr__(self):
        return '<User(email="%s")>' % self.email

    @staticmethod
    def create_user(
            email: str, password: str, crypted: bool, user_id: int = None,
            affiliation: Dict[str, str] = None, token: str = None, generate_token: bool = True,
            **kwargs):
        repo_db = infrastructure.repository_db
        repo_db.begin()
        try:
            if affiliation is not None:
                affiliation = Affiliation(**affiliation)
                repo_db.add(affiliation)

            user = User(email=email, user_id=user_id, affiliation=affiliation, **kwargs)
            repo_db.add(user)
            user.set_password(password, crypted)

            # TODO this has to change, e.g. trade for JWTs
            if token is None and generate_token:
                token = ''.join(random.choices(User._token_chars, k=64))
            if token is not None:
                repo_db.add(Session(token=token, user=user))

            repo_db.commit()
            return user
        except Exception as e:
            repo_db.rollback()
            utils.get_logger('__name__').error('could not create user', email=email, exc_info=e)
            raise e

    def update(self, crypted: bool = True, password: str = None, **kwargs):
        repo_db = infrastructure.repository_db
        repo_db.begin()
        try:
            if password is not None:
                self.set_password(password, crypted=crypted)

            for key in kwargs:
                setattr(self, key, kwargs.get(key))

            repo_db.commit()
        except Exception as e:
            repo_db.rollback()
            utils.get_logger('__name__').error(
                'could not edit user', email=self.email, user_id=self.user_id, exc_info=e)
            raise e

    def _verify_password(self, password):
        return bcrypt.verify(password, self.password)

    @staticmethod
    def from_user_id(user_id) -> 'User':
        return infrastructure.repository_db.query(User).get(user_id)

    def get_auth_token(self):
        repo_db = infrastructure.repository_db
        session = repo_db.query(Session).filter_by(user_id=self.user_id).first()

        if not session:
            # No session, user probably not logged in at NOMAD-coe repository GUI
            repo_db.begin()
            try:
                # TODO this has to change, e.g. trade for JWTs
                token = ''.join(random.choices(User._token_chars, k=64))
                session = Session(token=token, user=self)
                repo_db.add(session)

                repo_db.commit()
            except Exception as e:
                repo_db.rollback()
                utils.get_logger('__name__').error(
                    'could not generate token for user', email=self.email, user_id=self.user_id,
                    exc_info=e)
                raise e

        return session.token.encode('utf-8')

    def get_signature_token(self, expiration=10):
        """
        Genertes ver short term JWT token that can be used to sign download URLs.

        Returns: Tuple with token and expiration datetime
        """
        expires_at = datetime.datetime.utcnow() + datetime.timedelta(seconds=expiration)
        token = jwt.encode(
            dict(user=self.email, exp=expires_at),
            config.services.api_secret, 'HS256').decode('utf-8')
        return token, expires_at

    def set_password(self, password: str, crypted: bool):
        """
        Sets the users password. With ``crypted=True`` password is supposed to
        be already bcrypted and 2y-indented.
        """
        if password is None:
            return

        if crypted:
            self.password = password
        else:
            password_hash = bcrypt.encrypt(password, ident='2y')
            self.password = password_hash

    @property
    def token(self):
        return self.get_auth_token().decode('utf-8')

    @property
    def is_admin(self) -> bool:
        return self.email == 'admin'

    @staticmethod
    def verify_user_password(email, password):
        if email is None or password is None or email == '' or password == '':
            return None

        repo_db = infrastructure.repository_db
        user = repo_db.query(User).filter_by(email=email).first()
        if not user:
            return None

        if user._verify_password(password):
            return user
        else:
            raise LoginException('Wrong password')

    @staticmethod
    def verify_auth_token(token):
        if token is None or token == '':
            return None

        repo_db = infrastructure.repository_db
        session = repo_db.query(Session).filter_by(token=token).first()
        if session is None:
            return None

        user = repo_db.query(User).filter_by(user_id=session.user_id).first()
        assert user, 'User in sessions must exist.'
        return user

    @staticmethod
    def verify_signature_token(token):
        """
        Verifies the given JWT token. This should be used to verify URLs signed
        with a short term signature token (see :func:`get_signature_token`)
        """
        try:
            decoded = jwt.decode(token, config.services.api_secret, algorithms=['HS256'])
            repo_db = infrastructure.repository_db
            user = repo_db.query(User).filter_by(email=decoded['user']).first()
            if user is None:
                raise LoginException('Token signed for invalid user')

            return user
        except KeyError:
            raise LoginException('Token with invalid/unexpected payload')
        except jwt.ExpiredSignatureError:
            raise LoginException('Expired token')
        except jwt.InvalidTokenError:
            raise LoginException('Invalid token')

    def to_popo(self) -> utils.POPO:
        popo = utils.POPO(
            id=self.user_id,
            first_name=self.first_name,
            last_name=self.last_name,
            email=self.email)
        if self.affiliation is not None:
            popo.update(affiliation=dict(
                name=self.affiliation.name,
                address=self.affiliation.address))

        return popo


def ensure_test_user(email):
    """
    Allows tests to make sure that the default test users exist in the database.
    Returns:
        The user as :class:`User` instance.
    """
    repo_db = infrastructure.repository_db
    existing = repo_db.query(User).filter_by(email=email).first()
    assert existing, 'Test user %s does not exist.' % email

    session = repo_db.query(Session).filter_by(
        user_id=existing.user_id).first()
    assert session, 'Test user %s has no session.' % email
    assert session.token == existing.first_name.lower(), 'Test user %s session has unexpected token.' % email

    return existing


def admin_user():
    """
    Returns the admin user, a special user with `user_id==0`.
    Its password is part of :mod:`nomad.config`.
    """
    repo_db = infrastructure.repository_db
    admin = repo_db.query(User).filter_by(user_id=0).first()
    assert admin, 'Admin user does not exist.'
    return admin
