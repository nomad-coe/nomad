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

from typing_extensions import Annotated
from typing import Optional, Dict, Union, List, Literal, Any
from pydantic import BaseModel, Field, root_validator
import os.path
import pkgutil
import yaml

from .models import Options


class PluginBase(BaseModel):
    '''
    Base model for a NOMAD plugin.

    This should not be used. Plugins should instantiate concrete Plugin models like
    Parser or Schema.
    '''
    name: str = Field(description='A short descriptive human readable name for the plugin.')
    description: Optional[str] = Field(description='A human readable description of the plugin.')
    plugin_documentation_url: Optional[str] = Field(description='The URL to the plugins main documentation page.')
    plugin_source_code_url: Optional[str] = Field(description='The URL of the plugins main source code repository.')


class PythonPluginBase(PluginBase):
    '''
    A base model for NOMAD plugins that are implemented in Python.
    '''
    python_package: str = Field(description='''
        Name of the python package that contains the plugin code and a
        plugin metadata file called `nomad_plugin.yaml`.
    ''')

    @classmethod
    def _add_metadata(cls, metadata: Dict[str, Any], values: Dict[str, Any]):
        values.update(metadata)

    @root_validator(pre=True)
    def load_metadata(cls, values: Dict[str, Any]):  # pylint: disable=no-self-argument
        python_package = values.get('python_package')
        if not python_package:
            raise ValueError('Python plugins must provide a python_package.')

        try:
            # We manually look for the package to avoid the circlular imports that
            # all "official" methods (importlib, pkgutil) will cause. This will also
            # make the config import faster.
            # We try to deduce the package path form the top-level package
            package_path_segments = python_package.split('.')
            root_package = package_path_segments[0]
            package_dirs = package_path_segments[1:]
            package_path = os.path.join(
                os.path.dirname(pkgutil.get_loader(root_package).get_filename()),  # type: ignore
                *package_dirs)
            if not os.path.isdir(package_path):
                # We could not find it this way. Let's try to official way
                package_path = os.path.dirname(pkgutil.get_loader(python_package).get_filename())  # type: ignore
        except Exception as e:
            raise ValueError(f'The python package {python_package} cannot be loaded.', e)

        metadata_path = os.path.join(package_path, 'nomad_plugin.yaml')
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, 'r', encoding='UTF-8') as f:
                    metadata = yaml.load(f, Loader=yaml.FullLoader)
            except Exception as e:
                raise ValueError(f'Cannot load plugin metadata file {metadata_path}.', e)

            cls._add_metadata(metadata, values)
        return values


class Schema(PythonPluginBase):
    '''
    A Schema describes a NOMAD Python schema that can be loaded as a plugin.
    '''
    plugin_type: Literal['schema'] = Field('schema', description='''
        The type of the plugin. This has to be the string `schema` for schema plugins.
    ''')


class Parser(PythonPluginBase):
    '''
    A Parser describes a NOMAD parser that can be loaded as a plugin.

    The parser itself is references via `python_name`. For Parser instances `python_name`
    must refer to a Python class that has a `parse` function. The other properties are
    used to create a `MatchingParserInterface`. This comprises general metadata that
    allows users to understand what the parser is, and metadata used to decide if a
    given file "matches" the parser.
    '''

    # TODO the nomad_plugin.yaml for each parser needs some cleanup. The way parser metadata
    #      is presented in the UIs should be rewritten
    # TODO ideally we can somehow load parser plugin models lazily. Right now importing
    #      config will open all `nomad_plugin.yaml` files. But at least there is no python import
    #      happening.
    # TODO this should fully replace MatchingParserInterface
    # TODO most actual parser do not implement any abstract class. The Parser class has an
    #      abstract is_mainfile, which does not allow to separate parser implementation and plugin
    #      definition.

    plugin_type: Literal['parser'] = Field('parser', description='''
        The type of the plugin. This has to be the string `parser` for parser plugins.
    ''')

    parser_class_name: str = Field(description='''
        The fully qualified name of the Python class that implements the parser.
        This class must have a function `def parse(self, mainfile, archive, logger)`.
    ''')
    parser_as_interface: bool = Field(False, description='''
        By default the parser metadata from this config (and the loaded nomad_plugin.yaml)
        is used to instantiate a parser interface that is lazy loading the actual parser
        and performs the mainfile matching. If the parser interface matching
        based on parser metadata is not sufficient and you implemented your own
        is_mainfile parser method, this setting can be used to use the given
        parser class directly for parsing and matching.
    ''')

    mainfile_contents_re: Optional[str] = Field(description='''
        A regular expression that is applied the content of a potential mainfile.
        If this expression is given, the parser is only considered for a file, if the
        expression matches.
    ''')
    mainfile_name_re: str = Field(r'.*', description='''
        A regular expression that is applied the name of a potential mainfile.
        If this expression is given, the parser is only considered for a file, if the
        expression matches.
    ''')
    mainfile_mime_re: str = Field(r'text/.*', description='''
        A regular expression that is applied the mime type of a potential mainfile.
        If this expression is given, the parser is only considered for a file, if the
        expression matches.
    ''')
    mainfile_binary_header: Optional[bytes] = Field(description='''
        Matches a binary file if the given bytes are included in the file.
    ''')
    mainfile_binary_header_re: Optional[bytes] = Field(description='''
        Matches a binary file if the given binary regular expression bytes matches the
        file contents.
    ''')
    mainfile_alternative: bool = Field(False, description='''
        If True, the parser only matches a file, if no other file in the same directory
        matches a parser.
    ''')
    mainfile_contents_dict: Optional[dict] = Field(description='''
        Is used to match structured data files like JSON or HDF5.
    ''')
    supported_compressions: List[str] = Field([], description='''
        Files compressed with the given formats (e.g. xz, gz) are uncompressed and
        matched like normal files.
    ''')
    domain: str = Field('dft', description='''
        The domain value `dft` will apply all normalizers for atomistic codes. Deprecated.
    ''')
    level: int = Field(0, description='''
        The order by which the parser is executed with respect to other parsers.
    ''')

    code_name: Optional[str]
    code_homepage: Optional[str]
    code_category: Optional[str]
    metadata: Optional[dict] = Field(description='''
        Metadata passed to the UI. Deprecated. ''')

    def create_matching_parser_interface(self):
        if self.parser_as_interface:
            from nomad.parsing.parser import import_class
            Parser = import_class(self.parser_class_name)
            return Parser()

        from nomad.parsing.parser import MatchingParserInterface
        data = self.dict()
        del data['description']
        del data['python_package']
        del data['plugin_type']
        del data['parser_as_interface']
        del data['plugin_source_code_url']
        del data['plugin_documentation_url']

        return MatchingParserInterface(**data)


Plugin = Annotated[Union[Schema, Parser], Field(discriminator='plugin_type')]


class Plugins(Options):
    options: Dict[str, Plugin] = Field(dict(), description='The available plugin.')
