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

''' The beginning of a pylint plugin. Unfotunately it is kinda nonsensical without a
partnering mypy plugin. '''

import astroid
from astroid import MANAGER


annotation_names = {
    'MSection': {
        'a_test': '*'
    },
    'Section': {
        'a_elastic': 'nomad.metainfo.elastic_extension.ElasticDocument'
    },
    'Quantity': {
        'a_elastic': 'nomad.metainfo.elastic_extension.Elastic'
    }
}


def register(linter):
    # Needed for registering the plugin.
    pass


def transform(cls):
    for cls_name, annotations in annotation_names.items():
        if cls.name == cls_name:
            for name, type_spec in annotations.items():
                if type_spec == '*':
                    cls.locals[name] = [astroid.Instance()]
                else:
                    type_path = type_spec.split('.')
                    type_module = '.'.join(type_path[:-1])
                    type_name = type_path[-1]
                    module = MANAGER.ast_from_module_name(type_module)
                    cls.locals[name] = [cls.instantiate_class() for cls in module.lookup(type_name)[1]]


MANAGER.register_transform(astroid.ClassDef, transform)
