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

'''
API endpoint that deliver backend configuration details.
'''

from typing import Dict, Any, List, Optional
from datetime import datetime
from fastapi.routing import APIRouter
from pydantic.fields import Field
from pydantic.main import BaseModel

from nomad import config, normalizing
from nomad.utils import strip
from nomad.search import search
from nomad.parsing import parsers
from nomad.parsing.parsers import code_metadata
from nomad.app.v1.models import Aggregation, StatisticsAggregation
from nomad.metainfo.elasticsearch_extension import entry_type


router = APIRouter()
default_tag = 'info'


class MetainfoModel(BaseModel):
    all_package: str = Field(None, description=strip('''
        Name of the metainfo package that references all available packages, i.e.
        the complete metainfo.'''))

    root_section: str = Field(None, description=strip('''
        Name of the topmost section, e.g. section run for computational material science
        data.'''))


class StatisticsModel(BaseModel):
    n_entries: int = Field(None, description='Number of entries in NOMAD')
    n_uploads: int = Field(None, description='Number of uploads in NOMAD')
    n_quantities: int = Field(None, description='Accumulated number of quantities over all entries in the Archive')
    n_calculations: int = Field(None, description='Accumulated number of calculations, e.g. total energy calculations in the Archive')
    n_materials: int = Field(None, description='Number of materials in NOMAD')
    # TODO raw_file_size, archive_file_size


class CodeInfoModel(BaseModel):
    code_name: Optional[str] = Field(None, description='Name of the code or input format')
    code_homepage: Optional[str] = Field(None, description='Homepage of the code or input format')


class InfoModel(BaseModel):
    parsers: List[str]
    metainfo_packages: List[str]
    codes: List[CodeInfoModel]
    normalizers: List[str]
    statistics: StatisticsModel = Field(None, description='General NOMAD statistics')
    search_quantities: dict
    version: str
    deployment: str
    oasis: bool
    # TODO this should be removed in later releases, once most regular NOMAD users
    # should have switched to a new GUI version.
    git: dict = Field(None, description=strip('''
        A deprecated field that always contains an empty value to retain some compatibility
        with older GUIs.
    '''))


_statistics: Dict[str, Any] = None


def statistics():
    global _statistics
    if _statistics is None or datetime.now().timestamp() - _statistics.get('timestamp', 0) > 3600 * 24:
        _statistics = dict(timestamp=datetime.now().timestamp())
        search_response = search(aggregations=dict(statistics=Aggregation(statistics=StatisticsAggregation(
            metrics=['n_entries', 'n_materials', 'n_uploads', 'n_quantities', 'n_calculations']))))
        _statistics.update(**search_response.aggregations['statistics'].statistics.data)  # pylint: disable=no-member

    return _statistics


@router.get(
    '',
    tags=[default_tag],
    summary='Get information about the nomad backend and its configuration',
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    response_model=InfoModel)
async def get_info():
    ''' Return information about the nomad backend and its configuration. '''
    return InfoModel(**{
        'parsers': [
            key[key.index('/') + 1:]
            for key in parsers.parser_dict.keys()],
        'metainfo_packages': ['general', 'general.experimental', 'common', 'public'] + sorted([
            key[key.index('/') + 1:]
            for key in parsers.parser_dict.keys()]),
        'codes': [
            {'code_name': x['codeLabel'], 'code_homepage': x['codeUrl']}
            for x in sorted(code_metadata.values(), key=lambda info: info['codeLabel'].lower())
        ],
        'normalizers': [normalizer.__name__ for normalizer in normalizing.normalizers],
        'statistics': statistics(),
        'search_quantities': {
            s.qualified_name: {
                'name': s.qualified_name,
                'description': s.definition.description,
                'many': not s.definition.is_scalar
            }
            for s in entry_type.quantities.values()
            if 'optimade' not in s.qualified_name
        },
        'version': config.meta.version,
        'deployment': config.meta.deployment,
        'oasis': config.oasis.is_oasis,
        'git': {}
    })
