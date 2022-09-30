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

from nomad.metainfo import MTypes, Datetime, MEnum, Reference


validElnTypes = {
    'str': ['str'],
    'bool': ['bool'],
    'number': [x.__name__ for x in MTypes.num_python] + [f'np.{x.__name__}' for x in MTypes.num_numpy],  # type: ignore
    'datetime': ['Datetime'],
    'enum': ['{type_kind: Enum, type_data: [Operator, Responsible_person]}'],
    'user': ['User'],
    'author': ['Author'],
    'reference': ['']
}

validElnComponents = {
    'str': ['StringEditQuantity', 'FileEditQuantity', 'RichTextEditQuantity', 'EnumEditQuantity'],
    'bool': ['BoolEditQuantity'],
    'number': ['NumberEditQuantity', 'SliderEditQuantity'],
    'datetime': ['DateTimeEditQuantity'],
    'enum': ['EnumEditQuantity', 'AutocompleteEditQuantity', 'RadioEnumEditQuantity'],
    'user': ['AuthorEditQuantity'],
    'author': ['AuthorEditQuantity'],
    'reference': ['ReferenceEditQuantity']
}


def validate_eln_quantity_annotations(quantity):
    def assert_component(component_name, quantity_name, quantity_type, accepted_components):
        assert component_name in accepted_components, (
            f'The component {component_name} is not compatible with the quantity '
            f'{quantity_name} of the type {quantity_type}. '
            f'Accepted components: {", ".join(accepted_components)}.')

    if 'eln' not in quantity.m_annotations:
        return

    component = quantity.m_annotations['eln'].get('component', False)
    assert component, 'Quantity ELN annotation need to define a component'

    mtype = quantity.type
    name = quantity.name
    if isinstance(mtype, type):
        if mtype.__name__ == 'str':
            assert_component(component, name, mtype.__name__, validElnComponents['str'])
        elif mtype.__name__ == 'bool':
            assert_component(component, name, mtype.__name__, validElnComponents['bool'])
        elif mtype in MTypes.num_python:
            assert_component(component, name, mtype.__name__, validElnComponents['number'])
        elif mtype in MTypes.num_numpy:
            assert_component(component, name, f'np.{mtype.__name__}', validElnComponents['number'])
        elif mtype.__name__ == 'User':
            assert_component(component, name, mtype.__name__, validElnComponents['user'])
        elif mtype.__name__ == 'Author':
            assert_component(component, name, mtype.__name__, validElnComponents['author'])
    elif mtype == Datetime:
        assert_component(component, name, type(mtype).__name__, validElnComponents['datetime'])
    elif isinstance(mtype, MEnum):
        assert_component(component, name, type(mtype).__name__, validElnComponents['enum'])
    elif isinstance(mtype, Reference):
        target_class = mtype.target_section_def.section_cls
        if target_class.__name__ == 'User':
            assert_component(component, name, target_class.__name__, validElnComponents['user'])
        elif target_class.__name__ == 'Author':
            assert_component(component, name, target_class.__name__, validElnComponents['author'])
        else:
            assert_component(component, name, type(mtype).__name__, validElnComponents['reference'])
