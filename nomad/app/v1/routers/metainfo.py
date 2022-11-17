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
import copy
import datetime
from typing import Any, Dict

from fastapi import APIRouter, Path, status, HTTPException
from pydantic import BaseModel, Field

from nomad.app.v1.models import HTTPExceptionModel
from nomad.app.v1.utils import create_responses
from nomad.metainfo import Package
from nomad.metainfo.metainfo import MSection, Section, Quantity, Datetime, JSON
from nomad.metainfo.mongoengine_extension import MongoDocument, Mongo
from nomad.utils import strip, get_logger
from nomad import config

logger = get_logger(__name__)


class PackageDefinition(MSection):
    m_def = Section(a_mongo=MongoDocument())

    definition_id = Quantity(
        type=str,
        description='The sha1 based 40-digit long unique id for the package.',
        a_mongo=Mongo(primary_key=True, regex=r'^\w{40}$')
    )
    upload_id = Quantity(
        type=str,
        description='The upload id of the upload that contains the package.',
        a_mongo=Mongo()
    )
    qualified_name = Quantity(
        type=str,
        description='The qualified name for the package.',
        a_mongo=Mongo()
    )
    date_created = Quantity(
        type=Datetime,
        description='The date when the package is stored in the database.',
        a_mongo=Mongo()
    )
    package_definition = Quantity(
        type=JSON,
        description='Plain JSON representation (Python dict) of the package.',
        a_mongo=Mongo()
    )
    section_definition_ids = Quantity(
        type=str,
        shape=['*'],
        description='A list of section unique ids defined in this package.',
        a_mongo=Mongo(index=True, regex=r'^\w{40}$')
    )

    def __init__(self, package: Package, **kwargs):
        super().__init__()
        self.definition_id = package.definition_id
        if 'upload_id' in kwargs:
            self.upload_id = kwargs['upload_id']
            del kwargs['upload_id']
        self.qualified_name = package.qualified_name()
        self.date_created = datetime.datetime.utcnow()
        self.package_definition = package.m_to_dict(
            with_def_id=config.process.write_definition_id_to_archive, **kwargs)
        self.section_definition_ids = [section.definition_id for section in package.section_definitions]
        self.quantity_definition_ids = [
            quantity.definition_id for section in package.section_definitions for quantity in section.quantities]


def store_package_definition(package: Package, **kwargs):
    if package is None:
        return

    if PackageDefinition.m_def.a_mongo.objects(definition_id=package.definition_id).count() > 0:
        logger.info(f'Package {package.definition_id} already exists. Skipping.')
        return

    mongo_package = PackageDefinition(package, **kwargs)
    mongo_package.a_mongo.save()


#
# FastAPI router for the metainfo API.
#


router = APIRouter()

metainfo_tag = 'metainfo'

_bad_definition_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''Package not found. The given section definition is not contained in any packages.''')
}

_not_authorized_to_upload = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''Unauthorized. No credentials provided.''')
}


class PackageDefinitionResponse(BaseModel):
    section_definition_id: str = Field(None)
    data: Dict[str, Any] = Field(None)


def get_package_by_section_definition_id(section_definition_id: str) -> dict:
    packages = PackageDefinition.m_def.a_mongo.objects(section_definition_ids=section_definition_id)

    if packages.count() == 0:
        raise HTTPException(
            status.HTTP_404_NOT_FOUND,
            detail='Package not found. The given section definition is not contained in any packages.'
        )

    result = packages.first()

    pkg_definition = result.package_definition
    # add entry_id_based_name as a field which will be later used as the package name
    pkg_definition['entry_id_based_name'] = str(result.qualified_name)

    return copy.deepcopy(pkg_definition)


@router.get(
    '/{section_definition_id}', tags=[metainfo_tag],
    summary='Get the definition of package that contains the target id based section definition.',
    response_model=PackageDefinitionResponse,
    responses=create_responses(_bad_definition_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_package_definition(
        section_definition_id: str = Path(
            ..., regex=r'^\w{40}$',
            description='The section definition id to be used to retrieve package.')
):
    '''
    Retrieve the package that contains the target section.
    '''
    return {
        'section_definition_id': section_definition_id,
        'data': get_package_by_section_definition_id(section_definition_id)}
