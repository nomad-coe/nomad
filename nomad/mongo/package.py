#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import datetime
from re import fullmatch

from fastapi import HTTPException
from mongoengine import DateTimeField, DictField, Document, ListField, StringField
from pymongo.errors import DocumentTooLarge
from starlette import status

from nomad.metainfo import Package
from nomad.utils import get_logger

logger = get_logger(__name__)


class PackageDefinition(Document):
    id_pattern = r'^[a-f0-9]{40}$'

    snapshot_package_id = StringField(primary_key=True, regex=id_pattern, required=True)
    date_created = DateTimeField(default=datetime.datetime.now)
    entry_id = StringField(required=False)
    upload_id = StringField(required=False)
    qualified_name = StringField(required=True)
    package_definition = DictField(required=True)
    snapshot_section_ids = ListField(StringField(regex=id_pattern), default=None)

    meta = {'indexes': ['snapshot_section_ids']}

    @classmethod
    def has_package(cls, package_id: str):
        return cls.objects(snapshot_package_id=package_id).count() > 0

    @classmethod
    def has_section(cls, section_id: str):
        return cls.objects(snapshot_section_ids=section_id).count() > 0

    @classmethod
    def has_definition(cls, definition_id: str):
        if not fullmatch(PackageDefinition.id_pattern, definition_id):
            return False

        return cls.has_package(definition_id) or cls.has_section(definition_id)

    @classmethod
    def create_new(cls, package: Package, *, overwrite_existing: bool = True, **kwargs):
        if package is None or package.loaded_from_mongodb:
            return

        def _fields():
            return dict(
                entry_id=package.entry_id,
                upload_id=package.upload_id,
                qualified_name=package.qualified_name(),
                package_definition=package.m_to_dict(
                    **(dict(with_def_id=True) | kwargs)
                ),
                snapshot_section_ids=[
                    section.definition_id for section in package.section_definitions
                ],
                date_created=datetime.datetime.now(),
            )

        target = cls.objects(snapshot_package_id=package.definition_id)

        try:
            if target.count() == 0:
                cls(snapshot_package_id=package.definition_id, **_fields()).save()
            elif overwrite_existing:
                target.update_one(**{f'set__{k}': v for k, v in _fields().items()})
        except DocumentTooLarge:
            logger.warning(
                f'Package definition for package {package.qualified_name()} is too large to be stored in MongoDB.'
            )

    @classmethod
    def get_by(cls, snapshot_id: str):
        """
        Get the package definition that contains the given section definition ID.
        """
        for field in ('snapshot_section_ids', 'snapshot_package_id'):
            if package := cls.objects(**{field: snapshot_id}).first():
                result = package.to_mongo().to_dict()
                result['snapshot_package_id'] = result.pop('_id')
                return result

        raise HTTPException(
            status.HTTP_404_NOT_FOUND,
            detail='Package not found. The given section definition is not contained in any packages.',
        )
