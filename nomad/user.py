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

"""
Module with some prototypes/placeholder for future user management in nomad@FAIR.
It is currently based on the NOMAD-coe repository postgres API. This module allows
to authenticate users based on user password or session tokens. It allows to access
the user data like names and user_id.

This implementation is based on SQLAlchemy. There are model classes that represent
entries in the *users* and *session* tables.

.. autoclass:: User
    :members:
    :undoc-members:

.. autoclass:: Session
    :members:
    :undoc-members:

.. autofunction:: ensure_test_user
"""

from passlib.hash import bcrypt
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base

from nomad import infrastructure


Base = declarative_base()


class Session(Base):  # type: ignore
    __tablename__ = 'sessions'

    token = Column(String, primary_key=True)
    user_id = Column(String)


class LoginException(Exception):
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
    firstname = Column(String)
    lastname = Column(String)
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

    def get_auth_token(self):
        repository_db = infrastructure.repository_db
        session = repository_db.query(Session).filter_by(user_id=self.user_id).first()
        if not session:
            raise LoginException('No session, user probably not logged in at NOMAD-coe repository GUI')

        return session.token.encode('utf-8')

    @staticmethod
    def verify_user_password(email, password):
        repository_db = infrastructure.repository_db
        user = repository_db.query(User).filter_by(email=email).first()
        if not user:
            return None

        if user._verify_password(password):
            return user
        else:
            raise LoginException('Wrong password')

    @staticmethod
    def verify_auth_token(token):
        repository_db = infrastructure.repository_db
        session = repository_db.query(Session).filter_by(token=token).first()
        if session is None:
            return None

        user = repository_db.query(User).filter_by(user_id=session.user_id).first()
        assert user, 'User in sessions must exist.'
        return user


def ensure_test_user(email):
    """
    Allows tests to make sure that the default test users exist in the database.
    Returns:
        The user as :class:`User` instance.
    """
    existing = infrastructure.repository_db.query(User).filter_by(
        email=email).first()
    assert existing, 'Test user %s does not exist.' % email

    session = infrastructure.repository_db.query(Session).filter_by(
        user_id=existing.user_id).first()
    assert session, 'Test user %s has no session.' % email
    assert session.token == email, 'Test user %s session has unexpected token.' % email

    return existing
