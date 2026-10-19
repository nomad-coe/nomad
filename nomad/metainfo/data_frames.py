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

import inspect
import types
from collections.abc import Iterable
from typing import Union, cast

import numpy as np
import xarray as xr
from pydantic import BaseModel

from nomad.metainfo import MSection, Package, Quantity, Section, SubSection, constraint
from nomad.units import ureg

m_package = Package()


class Values(MSection):
    name = Quantity(type=str)
    values = None
    values_ref = Quantity(type='Values', shape=[])
    spanned_dimensions = Quantity(type=int, shape=['*'])
    original_shape = Quantity(type=int, shape=['*'])

    def get_values(self, reshape: bool = True) -> np.ndarray:
        if self.values_ref:
            return self.values_ref.m_resolved().get_values()
        values = self.values
        if not isinstance(self.values, np.ndarray | ureg.Quantity):
            values = np.array(self.values)
        if reshape:
            return cast(np.ndarray, values).reshape(self.original_shape)
        return values

    def __init__(self, *args, **kwargs):
        values_ref = None
        values: list = []
        if len(args) == 0:
            pass
        elif len(args) == 1 and isinstance(args[0], list | np.ndarray | ureg.Quantity):
            values = args[0]
        elif len(args) == 1 and isinstance(args[0], (Values)):
            values_ref = args[0]
            values = None
        else:
            values = args

        original_shape = kwargs.pop('original_shape', None)
        if isinstance(values, np.ndarray):
            values_shape = values.shape
            quantity_shape = self.m_def.all_quantities['values'].shape[:-1]
            if len(values_shape) < len(quantity_shape):
                raise ValueError(
                    f'The quantity shape, {quantity_shape}, does not meet the '
                    f'lower-bound set by the values shape, {values_shape}'
                )
            flat_shape = values_shape[: len(quantity_shape)] + (-1,)
            values = values.reshape(flat_shape)
            if original_shape is None:
                original_shape = values_shape
        elif isinstance(values, Iterable):
            original_shape = [len(values)]

        super().__init__(
            values=values,
            values_ref=values_ref,
            original_shape=original_shape,
            **kwargs,
        )

    def xarray_attrs(self) -> dict[str, str]:
        return dict(
            units=self.m_def.all_quantities['values'].unit,
            long_name=self.m_def.all_quantities['values'].label,
            description=self.m_def.all_quantities['values'].description,
            iri=self.m_def.all_quantities['values'].iri,
        )


def _get_default_names(iterable: Iterable[Values]) -> list[str]:
    names = []
    for values in iterable:
        counter = 0
        while True:
            counter += 1
            unique_name = f'{values.m_def.name}_{counter}'
            if unique_name not in names:
                names.append(unique_name)
                break
    return names


def _get_names(iterable: Iterable[Values]) -> list[str]:
    default_names = _get_default_names(iterable)
    return [
        values.name if values.name else default
        for values, default in zip(iterable, default_names)
    ]


def _get_values(
    iterable: Iterable[Values], values: Union[str, 'ValuesTemplate']
) -> Values:
    return_values = None
    if isinstance(values, str):
        default_names = _get_default_names(iterable)
        for v, default in zip(iterable, default_names):
            if v.name == values or default == values:
                if return_values is not None:
                    raise ValueError(f'Multiple values matching {values}')
                return_values = v
        return return_values
    for v in iterable:
        if v.m_def == values.section_def:
            if return_values is not None:
                raise ValueError(f'Multiple values matching {values}')
            return_values = v
    return return_values


class DataFrame(MSection):
    fields = SubSection(section='Values', repeats=True)
    variables = SubSection(section='Values', repeats=True)

    def get_field(self, field: Union[str, 'ValuesTemplate']) -> Values:
        return _get_values(self.fields, field)

    def get_variable(self, variable: Union[str, 'ValuesTemplate']) -> Values:
        return _get_values(self.variables, variable)

    @constraint(warning=False)
    def check_dimensions(self):
        # TODO constrains that validate the soundness of field and variable dimensions
        pass

    @constraint(warning=False)
    def check_mandatory_fields_and_variables(self):
        data_frame_annotation = self.m_def.m_get_annotation(DataFrameAnnotation)
        if data_frame_annotation is not None:
            for index, field in enumerate(data_frame_annotation.mandatory_fields):
                assert index < len(self.fields), f'Mandatory field {field} missing'
                assert self.fields[index].m_def == field.section_def, (
                    f'Field {field} missing'
                )

            for index, variable in enumerate(data_frame_annotation.mandatory_variables):
                assert index < len(self.variables), (
                    f'Mandatory field {variable} missing'
                )
                assert self.variables[index].m_def == variable.section_def, (
                    f'Field {variable} missing'
                )

    def to_xarray(self) -> xr.Dataset:
        shape = []
        dims = []
        coords = {}
        var: Values
        for var, name in zip(self.variables, _get_names(self.variables)):
            if var.spanned_dimensions is None or len(var.spanned_dimensions) == 0:
                coord_dims = [name]
                shape.append(len(var.values))
                dims.append(name)
            elif len(var.spanned_dimensions) == 1:
                dim = var.spanned_dimensions[0]
                if dim >= len(shape):
                    shape.append(len(var.values))
                    dims.append(f'm_dim_{dim}')
                coord_dims = [f'm_dim_{dim}']
            else:
                raise NotImplementedError('Only one spanned dimension supported')
            coords[name] = (
                coord_dims,
                var.values,
                var.xarray_attrs(),
            )
        data_vars = {}
        field: Values
        for field, name in zip(self.fields, _get_names(self.fields)):
            data_vars[name] = (
                dims,
                cast(np.ndarray, field.values).reshape(shape),
                field.xarray_attrs(),
            )
        return xr.Dataset(
            data_vars=data_vars,
            coords=coords,
            attrs=dict(
                description=self.m_def.description,
                long_name=self.m_def.label,
            ),
        )

    def to_pandas(self):
        return self.to_xarray().to_dataframe()


def _get_package():
    package = inspect.currentframe().f_back.f_back.f_globals.get('m_package', None)
    assert package is not None, (
        'PhysicalQuantities have to be defined within a python package with global '
        'Package m_package variable'
    )
    assert isinstance(m_package, Package), 'm_package has to be a Package instance'
    return package


class ValuesTemplate:
    """
    A generator for quantities of a certain template with type, shape, unit, name, description, iri, etc.
    """

    def __init__(self, **kwargs):
        self.quantity = Quantity(**kwargs)
        assert self.quantity.name is not None, (
            'Values templates must be explicitly named'
        )

        class ValuesTemplate(Values):
            m_def = Section(name=self.quantity.name)
            values = self(name='values', shape=self.quantity.shape + ['*'])

        _get_package().section_definitions.append(ValuesTemplate.m_def)
        self.section_def = ValuesTemplate.m_def
        self.create = ValuesTemplate
        self.section_cls = ValuesTemplate

    def __call__(self, **kwargs):
        # Make a deep copy of the quantity via m_from_dict(m_to_dict)
        quantity = Quantity.m_from_dict(self.quantity.m_to_dict())
        quantity.m_update(**kwargs)
        return quantity


class DataFrameAnnotation(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    mandatory_fields: list[ValuesTemplate]
    mandatory_variables: list[ValuesTemplate]

    def dict(self, *args, **kwargs):
        return self.model_dump(*args, **kwargs)

    def model_dump(self, *args, **kwargs):
        return dict(
            mandatory_fields=[
                field.section_def.qualified_name() for field in self.mandatory_fields
            ],
            mandatory_variables=[
                variable.section_def.qualified_name()
                for variable in self.mandatory_variables
            ],
        )


class DataFrameTemplate:
    """
    A generator for data frames with specific mandatory fields and default variables.
    """

    def __init__(
        self,
        mandatory_fields: list[ValuesTemplate],
        mandatory_variables: list[ValuesTemplate] = [],
        **kwargs,
    ):
        self.sub_section = SubSection(**kwargs)
        self.fields = mandatory_fields
        self.variables = mandatory_variables

        assert self.sub_section.name is not None, (
            'DataFrame templates must be explicitly named'
        )

        class DataFrameTemplate(DataFrame):
            m_def = Section(name=self.sub_section.name)

            # TODO validation that default fields and variables are actually present

        DataFrameTemplate.m_def.m_annotations['data_frame'] = DataFrameAnnotation(
            mandatory_fields=mandatory_fields,
            mandatory_variables=mandatory_variables,
        )

        _get_package().section_definitions.append(DataFrameTemplate.m_def)
        self.create = DataFrameTemplate
        self.section_cls = DataFrameTemplate
        self.section_def = DataFrameTemplate.m_def
        self.sub_section.section = self.section_def

    def __call__(self, **kwargs):
        sub_section = self.sub_section.m_copy()
        sub_section.m_update(**kwargs)

        def __init_metainfo__(self):
            # TODO here we can add a more specialised section def to the caller
            # definition (e.g. MySection) as an inner_section_definition
            pass

        sub_section.__init_metainfo__ = types.MethodType(__init_metainfo__, sub_section)
        return sub_section


m_package.__init_metainfo__()
