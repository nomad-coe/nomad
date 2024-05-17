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

import sys

from nomad.metainfo import Environment
from .data import User, Author, user_reference, author_reference
from .datamodel import (
    Dataset,
    EditableUserMetadata,
    AuthLevel,
    MongoUploadMetadata,
    MongoEntryMetadata,
    MongoSystemMetadata,
    EntryMetadata,
    EntryArchive,
)
from .optimade import OptimadeEntry, Species
from .metainfo import m_env
from .results import Results
from .data import EntryData, ArchiveSection, Schema
from nomad.config.models.plugins import SchemaPackageEntryPoint
from .context import Context, ClientContext, ServerContext

m_env.m_add_sub_section(
    Environment.packages, sys.modules['nomad.datamodel.datamodel'].m_package
)  # type: ignore
m_env.m_add_sub_section(
    Environment.packages, sys.modules['nomad.datamodel.optimade'].m_package
)  # type: ignore


_all_metainfo_environment = None


def all_metainfo_packages():
    """
    Returns an Environment with all available Python metainfo packages. This will
    import all plugins, if they are not already imported.
    """
    from nomad.metainfo import Package, Environment
    from nomad.datamodel import EntryArchive
    from nomad.config.models.plugins import PythonPluginBase

    # Due to lazyloading plugins, we need to explicitly
    # import plugin's python packages if we want to assure that their
    # metainfo package is loaded.
    from nomad.config import config

    config.load_plugins()
    for entry_point in config.plugins.entry_points.filtered_values():
        if isinstance(entry_point, PythonPluginBase):
            entry_point.import_python_package()
        if isinstance(entry_point, SchemaPackageEntryPoint):
            entry_point.load()

    # Importing the parsers will also make sure that related schemas will be imported
    # even if they are not part of the plugin's python package as this will import
    # the parser class and not just the package.
    # TODO in practice this might not be necessary.
    from nomad.parsing.parsers import import_all_parsers

    import_all_parsers()

    # Create the ES mapping to populate ES annotations with search keys.
    from nomad.search import entry_type

    if not entry_type.mapping:
        entry_type.create_mapping(EntryArchive.m_def)

    # TODO we call __init_metainfo__() for all packages where this has been forgotten
    # by the package author. Ideally this would not be necessary and we fix the
    # actual package definitions.
    for module_key in sorted(list(sys.modules)):
        pkg: Package = getattr(sys.modules[module_key], 'm_package', None)
        if pkg is not None and isinstance(pkg, Package):
            if pkg.name not in Package.registry:
                pkg.__init_metainfo__()

    global _all_metainfo_environment
    if not _all_metainfo_environment:
        _all_metainfo_environment = Environment()
        for package in Package.registry.values():
            _all_metainfo_environment.m_add_sub_section(Environment.packages, package)

    return _all_metainfo_environment
