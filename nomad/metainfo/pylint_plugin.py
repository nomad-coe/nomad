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

''' The beginning of a pylint plugin. Unfotunately it is kinda nonsensical without a
partnering mypy plugin. '''

import astroid
from astroid import MANAGER


annotation_names = {
    'MSection': {
        'a_test': '*'
    },
    'Section': {
        'a_elastic': 'nomad.metainfo.elastic_extension.ElasticDocument',
        'a_mongo': 'nomad.metainfo.mongoengine_extension.MongoDocument'
    },
    'Quantity': {
        'a_elastic': 'nomad.metainfo.elastic_extension.Elastic'
    }
}


def register(linter):
    # Needed for registering the plugin.
    pass


def transform(cls):
    ''' Transforms annotation fields for known annotation classes. '''
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


def is_derived(func):
    ''' Check if this is call to the derived decorator. '''
    decorators = func.decorators.nodes if func.decorators else []
    for decorator_node in decorators:
        if isinstance(decorator_node, astroid.Call):
            if decorator_node.func and isinstance(decorator_node.func, astroid.Name):
                return decorator_node.func.name == 'derived'

    return False


def derived_transform(node, context=None):
    '''
    The derived decorator produces a Quantity. Pylint does not infer this on its own.
    We change the inferred type of a @derived call to a Quantity instance here.
    '''
    module = MANAGER.ast_from_module_name('nomad.metainfo.metainfo')
    class_defs = [cls.instantiate_class() for cls in module.lookup('Quantity')[1]]
    return iter([class_defs[0].instantiate_class()])


MANAGER.register_transform(astroid.ClassDef, transform)
MANAGER.register_transform(astroid.FunctionDef, astroid.inference_tip(derived_transform), is_derived)
