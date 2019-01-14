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

.. autoclass:: UploadMetaData
    :members:
.. autoclass:: Upload
    :members:
    :undoc-members:
.. autoclass:: Calc
    :members:
    :undoc-members:
"""

from typing import Type
import datetime
from sqlalchemy import Column, Integer, String, Boolean, DateTime, ForeignKey
from sqlalchemy.orm import relationship

from nomad import utils, infrastructure, datamodel, files

from .user import User
from .calc import Calc
from .base import Base, CalcMetaData, UserMetaData, StructRatio, CodeVersion, Spacegroup, \
    CalcSet, Citation


topic_code = 220
topic_atoms = 10
topic_system_type = 50
topic_xc_treatment = 75
topic_crystal_system = 90
topic_basis_set_type = 80


class UploadMetaData:
    """
    Utility class that provides per upload meta data and overriding per calculation
    meta data. For a given *mainfile* data is first read from the `calculations` key
    (a list of calculation dict with a matching `mainfile` key), before it is read
    from :param:`meta_data_dict` it self.

    The class is used to deal with user provided meta-data.

    Arguments:
        meta_data_dict: The python dict with the meta-data.
    """
    def __init__(self, meta_data_dict: dict) -> None:
        self._upload_data = meta_data_dict
        self._calc_data: dict = {
            calc['mainfile']: calc
            for calc in self._upload_data.get('calculations', [])}

    def get(self, mainfile: str) -> dict:
        return self._calc_data.get(mainfile, self._upload_data)


class Upload(Base, datamodel.Upload):  # type: ignore
    __tablename__ = 'uploads'

    coe_upload_id = Column('upload_id', Integer, primary_key=True, autoincrement=True)
    upload_name = Column(String)
    user_id = Column(Integer, ForeignKey('users.user_id'))
    is_processed = Column(Boolean)
    created = Column(DateTime)

    user = relationship('User')
    calcs = relationship('Calc')

    @classmethod
    def load_from(cls, obj):
        return Upload.from_upload_id(str(obj.upload_id))

    @staticmethod
    def from_upload_id(upload_id: str) -> 'Upload':
        repo_db = infrastructure.repository_db
        uploads = repo_db.query(Upload).filter_by(upload_name=upload_id)
        assert uploads.count() <= 1, 'Upload id/name must be unique'
        return uploads.first()

    @property
    def upload_id(self) -> str:
        return self.upload_name

    @property
    def uploader(self) -> 'User':
        return self.user

    @property
    def upload_time(self) -> Type[datetime.datetime]:
        return self.created

    @staticmethod
    def add(upload: datamodel.Upload, meta_data: dict = {}) -> int:
        """
        Add the upload to the NOMAD-coe repository db. It creates an
        uploads-entry, respective calculation and property entries. Everything in one
        transaction.

        Triggers and updates the NOMAD-coe repository elastic search index after
        success (TODO).

        Arguments:
            upload: The upload to add.
            upload_meta_data: A dictionary with additional meta data (e.g. user provided
                meta data) that should be added to upload and calculations.
        """
        upload_meta_data = UploadMetaData(meta_data)
        repo_db = infrastructure.repository_db
        repo_db.begin()

        logger = utils.get_logger(__name__, upload_id=upload.upload_id)

        result = None

        try:
            # create upload
            coe_upload = Upload(
                upload_name=upload.upload_id,
                created=meta_data.get('_upload_time', upload.upload_time),
                user=upload.uploader,
                is_processed=True)
            repo_db.add(coe_upload)

            # add calculations and metadata
            has_calcs = False
            for calc in upload.calcs:
                has_calcs = True
                coe_upload._add_calculation(calc.to(files.Calc), upload_meta_data.get(calc.mainfile))

            # commit
            if has_calcs:
                # empty upload case
                repo_db.commit()
                result = coe_upload.coe_upload_id
            else:
                repo_db.rollback()
        except Exception as e:
            logger.error('Unexpected exception.', exc_info=e)
            repo_db.rollback()
            raise e

        # TODO trigger index update
        pass

        return result

    def _add_calculation(self, calc: files.Calc, calc_meta_data: dict) -> None:
        repo_db = infrastructure.repository_db

        # table based properties
        coe_calc = Calc(
            coe_calc_id=calc_meta_data.get('_pid', None),
            checksum=calc_meta_data.get('_checksum', calc.calc_id),
            upload=self)
        repo_db.add(coe_calc)

        program_version = calc.program_version  # TODO shorten version names
        code_version = repo_db.query(CodeVersion).filter_by(content=program_version).first()
        if code_version is None:
            code_version = CodeVersion(content=program_version)
            repo_db.add(code_version)

        metadata = CalcMetaData(
            calc=coe_calc,
            added=calc_meta_data.get('_upload_time', self.upload_time),
            chemical_formula=calc.chemical_composition,
            filenames=('[%s]' % ','.join(['"%s"' % filename for filename in calc.files])).encode('utf-8'),
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
            permission=(1 if calc_meta_data.get('with_embargo', False) else 0))
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
        owner_user_id = calc_meta_data.get('_uploader', int(self.user_id))
        coe_calc.owners.append(repo_db.query(User).get(owner_user_id))

        for coauthor_id in calc_meta_data.get('coauthors', []):
            coe_calc.coauthors.append(repo_db.query(User).get(coauthor_id))

        for shared_with_id in calc_meta_data.get('shared_with', []):
            coe_calc.shared_with.append(repo_db.query(User).get(shared_with_id))

        # datasets
        for dataset_id in calc_meta_data.get('datasets', []):
            dataset = CalcSet(parent_calc_id=dataset_id, children_calc_id=coe_calc.coe_calc_id)
            repo_db.add(dataset)

        # references
        for reference in calc_meta_data.get('references', []):
            citation = repo_db.query(Citation).filter_by(
                value=reference,
                kind='EXTERNAL').first()

            if citation is None:
                citation = Citation(value=reference, kind='EXTERNAL')
                repo_db.add(citation)

            coe_calc.citations.append(citation)
