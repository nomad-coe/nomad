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
Interface to the NOMAD-coe repository postgres database. This implementation is based on
SQLAlchemy. There are model classes that represent entries in the *users* and *session*
tables.

This module allows to authenticate users based on user password or session tokens.
It allows to access the user data like names and user_id.

.. autoclass:: User
    :members:
    :undoc-members:

.. autoclass:: Session
    :members:
    :undoc-members:

.. autofunction:: ensure_test_user

This module also provides functionality to add parsed calculation data to the db:

.. autofunction:: add_upload
"""

import itertools
from passlib.hash import bcrypt
from sqlalchemy import Column, Integer, String, Boolean, DateTime, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import BYTEA

from nomad import utils, infrastructure
from nomad.repo import RepoCalc


Base = declarative_base()


def add_upload(upload, restricted: bool) -> int:
    """
    Add the processed upload to the NOMAD-coe repository db. It creates an
    uploads-entry, respective calculation and property entries. Everything in one
    transaction. Triggers an updates the NOMAD-coe repository elastic search index after
    success.

    TODO deal with the restricted parameter
    """
    repo_db = infrastructure.repository_db
    repo_db.begin()

    logger = utils.get_logger(
        __name__,
        upload_id=upload.upload_id,
        upload_hash=upload.upload_hash)

    result = None

    try:
        # create upload
        coe_upload = Upload(
            upload_name=upload.upload_hash,
            created=upload.upload_time,
            user_id=int(upload.user_id),
            is_processed=True)
        repo_db.add(coe_upload)

        # add calculations and metadata
        has_calcs = False
        for repo_calc in RepoCalc.upload_calcs(upload.upload_id):
            has_calcs = True
            add_calculation(upload, coe_upload, repo_calc, restricted)

        # commit
        if has_calcs:
            # empty upload case
            repo_db.commit()
            result = coe_upload.upload_id
        else:
            repo_db.rollback()
    except Exception as e:
        logger.error('Unexpected exception.', exc_info=e)
        repo_db.rollback()
        raise e

    # trigger index update
    pass

    return result


def add_calculation(upload, coe_upload, calc: RepoCalc, restricted: bool) -> None:
    repo_db = infrastructure.repository_db

    # table based properties
    coe_calc = Calc(checksum=calc.calc_hash, upload=coe_upload)
    repo_db.add(coe_calc)

    program_version = calc.program_version  # TODO shorten version names
    code_version = repo_db.query(CodeVersion).filter_by(content=program_version).first()
    if code_version is None:
        code_version = CodeVersion(content=program_version)
        repo_db.add(code_version)

    filenames = itertools.chain([calc.mainfile], calc.aux_files)

    metadata = CalcMetaData(
        calc=coe_calc,
        added=upload.upload_time,
        chemical_formula=calc.chemical_composition,
        filenames=('[%s]' % ','.join(['"%s"' % filename for filename in filenames])).encode('utf-8'),
        location=calc.mainfile,
        version=code_version)
    repo_db.add(metadata)

    struct_ratio = StructRatio(
        calc=coe_calc,
        chemical_formula=calc.chemical_composition,
        formula_units=1, nelem=1)
    repo_db.add(struct_ratio)

    user_metadata = UserMetaData(
        calc=coe_calc,
        permission=0 if not restricted else 1)
    repo_db.add(user_metadata)

    spacegroup = Spacegroup(
        calc=coe_calc,
        n=int(calc.space_group_number)
    )
    repo_db.add(spacegroup)

    # topic based properties
    coe_calc.set_value(topic_code, calc.program_name)
    for atom in set(calc.atom_species):
        coe_calc.set_value(topic_atoms, str(atom))  # TODO atom label not number
    coe_calc.set_value(topic_system_type, calc.system_type)
    coe_calc.set_value(topic_xc_treatment, calc.XC_functional_name)  # TODO function->treatment
    coe_calc.set_value(topic_crystal_system, calc.crystal_system)
    coe_calc.set_value(topic_basis_set_type, calc.basis_set_type)


class Calc(Base):  # type: ignore
    __tablename__ = 'calculations'

    calc_id = Column(Integer, primary_key=True, autoincrement=True)
    origin_id = Column(Integer, ForeignKey('uploads.upload_id'))
    upload = relationship('Upload')
    checksum = Column(String)

    def set_value(self, topic_cid: int, value: str) -> None:
        if value is None:
            return

        repo_db = infrastructure.repository_db
        topic = repo_db.query(Topics).filter_by(topic=value).first()
        if not topic:
            topic = Topics(cid=topic_cid, topic=value)
            repo_db.add(topic)

        tag = Tag(calc=self, topic=topic)
        repo_db.add(tag)


class CalcMetaData(Base):  # type: ignore
    __tablename__ = 'metadata'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    added = Column(DateTime)
    chemical_formula = Column(String)
    filenames = Column(BYTEA)
    location = Column(String)
    version_id = Column(Integer, ForeignKey('codeversions.version_id'))
    version = relationship('CodeVersion')


class UserMetaData(Base):  # type: ignore
    __tablename__ = 'user_metadata'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    permission = Column(Integer)


class StructRatio(Base):  # type: ignore
    __tablename__ = 'struct_ratios'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    formula_units = Column(Integer)
    nelem = Column(Integer)
    chemical_formula = Column(String)


class CodeVersion(Base):  # type: ignore
    __tablename__ = 'codeversions'

    version_id = Column(Integer, primary_key=True, autoincrement=True)
    content = Column(String)


class Spacegroup(Base):  # type: ignore
    __tablename__ = 'spacegroups'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    n = Column(Integer)


class Tag(Base):  # type: ignore
    __tablename__ = 'tags'
    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    tid = Column(Integer, ForeignKey('topics.tid'), primary_key=True)
    topic = relationship('Topics')

    def __repr__(self):
        return '<Tag(calc_id="%d", tid="%d)>' % (self.calc_id, self.tid)


topic_code = 220
topic_atoms = 10
topic_system_type = 50
topic_xc_treatment = 75
topic_crystal_system = 90
topic_basis_set_type = 80


class Topics(Base):  # type: ignore
    __tablename__ = 'topics'
    tid = Column(Integer, primary_key=True, autoincrement=True)
    cid = Column(Integer)
    topic = Column(String)


class Upload(Base):  # type: ignore
    __tablename__ = 'uploads'

    upload_id = Column(Integer, primary_key=True, autoincrement=True)
    upload_name = Column(String)
    user_id = Column(Integer, ForeignKey('users.user_id'))
    user = relationship('User')
    is_processed = Column(Boolean)
    created = Column(DateTime)


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
        repo_db = infrastructure.repository_db
        session = repo_db.query(Session).filter_by(user_id=self.user_id).first()
        if not session:
            raise LoginException('No session, user probably not logged in at NOMAD-coe repository GUI')

        return session.token.encode('utf-8')

    @staticmethod
    def verify_user_password(email, password):
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
        repo_db = infrastructure.repository_db
        session = repo_db.query(Session).filter_by(token=token).first()
        if session is None:
            return None

        user = repo_db.query(User).filter_by(user_id=session.user_id).first()
        assert user, 'User in sessions must exist.'
        return user


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
    assert session.token == email, 'Test user %s session has unexpected token.' % email

    return existing
