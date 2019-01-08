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
import enum
from passlib.hash import bcrypt
from sqlalchemy import Column, Integer, String, Boolean, DateTime, ForeignKey, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import BYTEA

from nomad import utils, infrastructure
from nomad.repo import RepoCalc


Base = declarative_base()


class UploadMetaData:
    def __init__(self, meta_data_dict):
        self._upload_data = meta_data_dict
        self._calc_data = {
            calc['mainfile']: calc
            for calc in self._upload_data.get('calculations', [])}

    def get(self, mainfile):
        return self._calc_data.get(mainfile, self._upload_data)


def add_upload(upload, meta_data: dict = {}) -> int:
    """
    Add the processed upload to the NOMAD-coe repository db. It creates an
    uploads-entry, respective calculation and property entries. Everything in one
    transaction.

    Triggers and updates the NOMAD-coe repository elastic search index after
    success (TODO).

    Arguments:
        upload: The upload to add
        upload_meta_data: A dictionary with additional meta data (e.g. user meta data)
            that should be added to upload and calculations.

    TODO deal with meta_data
    """
    upload_meta_data = UploadMetaData(meta_data)
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
            created=meta_data.get('_upload_time', upload.upload_time),
            user_id=int(upload.user_id),
            is_processed=True)
        repo_db.add(coe_upload)

        # add calculations and metadata
        has_calcs = False
        for repo_calc in RepoCalc.upload_calcs(upload.upload_id):
            has_calcs = True
            add_calculation(
                upload, coe_upload, repo_calc, upload_meta_data.get(repo_calc.mainfile))

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


def add_calculation(upload, coe_upload, calc: RepoCalc, calc_meta_data: dict) -> None:
    repo_db = infrastructure.repository_db

    # table based properties
    coe_calc = Calc(
        calc_id=calc_meta_data.get('_pid', None),
        checksum=calc_meta_data.get('_checksum', calc.calc_hash),
        upload=coe_upload)
    repo_db.add(coe_calc)

    program_version = calc.program_version  # TODO shorten version names
    code_version = repo_db.query(CodeVersion).filter_by(content=program_version).first()
    if code_version is None:
        code_version = CodeVersion(content=program_version)
        repo_db.add(code_version)

    filenames = itertools.chain([calc.mainfile], calc.aux_files)

    metadata = CalcMetaData(
        calc=coe_calc,
        added=calc_meta_data.get('_upload_time', upload.upload_time),
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
        label=calc_meta_data.get('comment', None),
        permission=1 if calc_meta_data.get('restricted', False) else 0)
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

    # user relations
    owner_user_id = calc_meta_data.get('uploader', int(upload.user_id))
    ownership = Ownership(calc_id=coe_calc.calc_id, user_id=owner_user_id)
    repo_db.add(ownership)

    for coauthor_id in calc_meta_data.get('coauthors', []):
        coauthorship = CoAuthorship(calc_id=coe_calc.calc_id, user_id=int(coauthor_id))
        repo_db.add(coauthorship)

    for shared_with_id in calc_meta_data.get('share_with', []):
        shareship = Shareship(calc_id=coe_calc.calc_id, user_id=int(shared_with_id))
        repo_db.add(shareship)

    # datasets
    for dataset_id in calc_meta_data.get('_datasets', []):
        dataset = CalcSet(parent_calc_id=dataset_id, children_calc_id=coe_calc.calc_id)
        repo_db.add(dataset)

    # references
    for reference in calc_meta_data.get('references', []):
        citation = repo_db.query(Citation).filter_by(
            label=reference,
            kind=CitationKind.EXTERNAL).first()

        if citation is None:
            citation = Citation(label=reference, kind=CitationKind.EXTERNAL)
            repo_db.add(citation)

        metadata_citation = MetaDataCitation(
            calc_id=coe_calc.calc_id,
            citation_id=citation.citation_id)

        repo_db.add(metadata_citation)


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
    label = Column(String)
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
        return '<Tag(calc_id="%d", tid="%d)>' % (int(self.calc_id), int(self.tid))


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


class Ownership(Base):  # type: ignore
    __tablename__ = 'ownerships'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    user_id = Column(Integer, ForeignKey('users.user_id'), primary_key=True)


class CoAuthorship(Base):  # type: ignore
    __tablename__ = 'coauthorships'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    user_id = Column(Integer, ForeignKey('users.user_id'), primary_key=True)


class Shareship(Base):  # type: ignore
    __tablename__ = 'shareships'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    user_id = Column(Integer, ForeignKey('users.user_id'), primary_key=True)


class CalcSet(Base):  # type: ignore
    __tablename__ = 'calcsets'

    parent_calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    children_calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)


class MetaDataCitation(Base):  # type: ignore
    __tablename__ = 'metadata_citations'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    citation_id = Column(Integer, ForeignKey('citations.citation_id'), primary_key=True)


class CitationKind(enum.Enum):
    INTERNAL = 1
    EXTERNAL = 2


class Citation(Base):  # type: ignore
    __tablename__ = 'citations'

    citation_id = Column(Integer, primary_key=True)
    value = Column(String)
    kind = Column(Enum(CitationKind))


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


def admin_user():
    repo_db = infrastructure.repository_db
    admin = repo_db.query(User).filter_by(user_id=1).first()
    assert admin, 'Admin user does not exist.'
    return admin
