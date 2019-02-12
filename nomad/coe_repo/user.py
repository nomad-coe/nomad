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

from passlib.hash import bcrypt
from sqlalchemy import Column, Integer, String
import datetime
import jwt

from nomad import infrastructure, config, utils

from .base import Base


class Session(Base):  # type: ignore
    __tablename__ = 'sessions'

    token = Column(String, primary_key=True)
    user_id = Column(String)


class LoginException(Exception):
    """ Exception that is raised if the user could not be logged in despite present
    credentials. """
    pass


class User(Base):  # type: ignore
    """
    SQLAlchemy model class that represents NOMAD-coe repository postgresdb *users*.
    Provides functions for authenticating via password or session token.

    It is not intended to create or update users. This should be done via the
    NOMAD-coe repository GUI.
    """
    __tablename__ = 'users'

    user_id = Column(Integer, primary_key=True)
    email = Column(String)
    first_name = Column(String, name='firstname')
    last_name = Column(String, name='lastname')
    affiliation = Column(String)
    password = Column(String)

    def __repr__(self):
        return '<User(email="%s")>' % self.email

    def _hash_password(self, password):
        assert False, 'Login functions are done by the NOMAD-coe repository GUI'
        # password_hash = bcrypt.encrypt(password, ident='2y')
        # self.password = password_hash

    def _verify_password(self, password):
        return bcrypt.verify(password, self.password)

    def _generate_auth_token(self, expiration=600):
        assert False, 'Login functions are done by the NOMAD-coe repository GUI'

    @staticmethod
    def from_user_id(user_id) -> 'User':
        return infrastructure.repository_db.query(User).get(user_id)

    def get_auth_token(self):
        repo_db = infrastructure.repository_db
        session = repo_db.query(Session).filter_by(user_id=self.user_id).first()
        if not session:
            raise LoginException('No session, user probably not logged in at NOMAD-coe repository GUI')

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
            else:
                return user
        except KeyError:
            raise LoginException('Token with invalid/unexpected payload')
        except jwt.ExpiredSignatureError:
            raise LoginException('Expired token')
        except jwt.InvalidTokenError:
            raise LoginException('Invalid token')

    def to_popo(self) -> utils.POPO:
        return utils.POPO(
            id=self.user_id,
            first_name=self.first_name,
            last_name=self.last_name,
            email=self.email,
            affiliation=self.affiliation)


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
    admin = repo_db.query(User).filter_by(user_id=1).first()
    assert admin, 'Admin user does not exist.'
    return admin
