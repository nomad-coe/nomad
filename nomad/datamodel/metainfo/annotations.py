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

from typing import List, Any, Union, Dict, Optional
from enum import Enum
from pydantic import Field, validator
import re

from pydantic.main import BaseModel

from nomad.utils import strip
from nomad.metainfo import AnnotationModel, MEnum, MTypes, Datetime, Reference, Quantity


class ELNComponentEnum(str, Enum):
    StringEditQuantity = 'StringEditQuantity'
    URLEditQuantity = 'URLEditQuantity'
    EnumEditQuantity = 'EnumEditQuantity'
    RadioEnumEditQuantity = 'RadioEnumEditQuantity'
    AutocompleteEditQuantity = 'AutocompleteEditQuantity'
    FileEditQuantity = 'FileEditQuantity'
    BoolEditQuantity = 'BoolEditQuantity'
    NumberEditQuantity = 'NumberEditQuantity'
    SliderEditQuantity = 'SliderEditQuantity'
    DateTimeEditQuantity = 'DateTimeEditQuantity'
    RichTextEditQuantity = 'RichTextEditQuantity'
    ReferenceEditQuantity = 'ReferenceEditQuantity'
    UserEditQuantity = 'UserEditQuantity'
    AuthorEditQuantity = 'AuthorEditQuantity'


valid_eln_types = {
    'str': ['str'],
    'bool': ['bool'],
    'number': [x.__name__ for x in MTypes.num_python] + [f'np.{x.__name__}' for x in MTypes.num_numpy],  # type: ignore
    'datetime': ['Datetime'],
    'enum': ['{type_kind: Enum, type_data: [Operator, Responsible_person]}'],
    'user': ['User'],
    'author': ['Author'],
    'reference': ['']
}


valid_eln_components = {
    'str': [
        ELNComponentEnum.StringEditQuantity,
        ELNComponentEnum.URLEditQuantity,
        ELNComponentEnum.FileEditQuantity,
        ELNComponentEnum.RichTextEditQuantity,
        ELNComponentEnum.EnumEditQuantity],
    'bool': [
        ELNComponentEnum.BoolEditQuantity],
    'number': [
        ELNComponentEnum.NumberEditQuantity,
        ELNComponentEnum.SliderEditQuantity],
    'datetime': [
        ELNComponentEnum.DateTimeEditQuantity],
    'enum': [
        ELNComponentEnum.EnumEditQuantity,
        ELNComponentEnum.AutocompleteEditQuantity,
        ELNComponentEnum.RadioEnumEditQuantity],
    'user': [
        ELNComponentEnum.AuthorEditQuantity],
    'author': [
        ELNComponentEnum.AuthorEditQuantity],
    'reference': [
        ELNComponentEnum.ReferenceEditQuantity]
}


class Filter(BaseModel):
    ''' A filter defined by an include list or and exclude list of the quantities or subsections. '''

    include: Optional[List[str]] = Field(
        None, description=strip('''
            The list of quantity or subsection names to be included.
        '''))
    exclude: Optional[List[str]] = Field(
        None, description=strip('''
            The list of quantity or subsection names to be excluded.
        '''))


class SectionProperties(BaseModel):
    ''' A filter defined by an include list or and exclude list of the quantities and subsections. '''

    visible: Optional[Filter] = Field(
        1, description=strip('''
            Defines the visible quantities and subsections.
        '''))
    editable: Optional[Filter] = Field(
        None, description=strip('''
            Defines the editable quantities and subsections.
        '''))
    order: Optional[List[str]] = Field(
        None, description=strip('''
            To customize the order of the quantities and subsections.
        '''))


class ELNAnnotation(AnnotationModel):
    '''
    These annotations control how data can be entered and edited.
    Use the key `eln` to add this annotations. For example:

    ```python
    class Sample(EntryData):
        sample_id = Quantity(type=str, a_eln=dict(component='StringEditQuantity'))`)
    ```

    or in YAML schemas:
    ```yaml
    Sample:
      quantities:
        sample_id:
          type: str
          m_annotations:
            eln:
              component: StringEditQuantity
    ```

    An `eln` annotation can be added to *section* and *quantity* definitions to different
    effects. In both cases, it controls how sections and quantities are represented in the GUI
    with different parameters; see below.

    The UI gives an overview about all ELN edit annotations and components
    [here](https://nomad-lab.eu/prod/v1/staging/gui/dev/editQuantity).
    '''

    component: ELNComponentEnum = Field(None, description='''
        The form field component that is used to make the annotated quantity editable.
        If no component is given, the quantity won't be editable. This can be used on quantities only.

        The supported values are:

        `StringEditQuantity`: For editing simple short string values.<br/>
        `URLEditQuantity`: For editing strings that are validated to be URLs.<br/>
        `EnumEditQuantity`: For Editing enum values. Uses a dropdown list with enum values. This component may be used for short enumerates.<br/>
        `RadioEnumEditQuantity`: For Editing enum values. Uses radio buttons.<br/>
        `AutocompleteEditQuantity`: For editing enum values. Uses an autocomplete form with dropdown list. This component may be used for longer enumerates.<br/>
        `FileEditQuantity`: For editing a reference to a file. Will allow to choose a file or upload a file.<br/>
        `BoolEditQuantity`: For editing boolean choices.<br/>
        `NumberEditQuantity`: For editing numbers with our without unit.<br/>
        `SliderEditQuantity`: For editing numbers with a horizontal slider widget.<br/>
        `DateTimeEditQuantity`: For editing datetimes.<br/>
        `RichTextEditQuantity`: For editing long styled text with a rich text editor.<br/>
        `ReferenceEditQuantity`: For editing references to other sections.<br/>
        `UserEditQuantity`: For entering user information. Lets you choose a nomad user or enter information manually.<br/>
        `AuthorEditQuantity`: For entering author information manually.
    ''')

    label: str = Field(None, description='Custom label for the quantity shown on the form field.')

    props: Dict[str, Any] = Field(None, description='''
        A dictionary with additional props that are passed to the  editcomponent.
    ''')

    default: Any = Field(None, description='''
        Prefills any set form field component with the given value. This is different
        from the quantities `default` property. The quantities default is not stored
        in the data; the default value is assumed if no other value is given. The
        ELN form field default value will be stored, even if not changed.
    ''')
    defaultDisplayUnit: str = Field(None, description='''
        Allows to define a default unit to initialize a `NumberEditQuantity` with. The
        unit has to be compatible with the unit of the annotation quantity and the annotated
        quantity must have a unit. Only applies to quantities and with
        `component=NumberEditQuantity`.
    ''')

    minValue: Union[int, float] = Field(None, description='''
        Allows to specify a minimum value for quantity annotations with number type.
        Will show an error, if outside numbers are entered. Only works on quantities and
        in conjunction with `component=NumberEditQuantity`.
    ''')
    maxValue: Union[int, float] = Field(None, description='''
        Allows to specify a maximum value for quantity annotations with number type.
        Will show an error, if outside numbers are entered. Only works on quantities and
        in conjunction with `component=NumberEditQuantity`.
    ''')

    hide: List[str] = Field(None, description='''
        The annotation "hide" is deprecated. Use "visible" key of "properties" annotation instead.
        Allows you to hide certain quantities from a section editor. Give a list
        of quantity names. Quantities must exist in the section that this annotation
        is added to. Can only be used in section annotations.
    ''', deprecated=True)

    overview: bool = Field(None, description='''
        Shows the annotation section on the entry's overview page. Can only be used on
        section annotations.''')

    lane_width: Union[str, int] = Field(None, description='''
        Value to overwrite the css width of the lane used to render the annotation
        section and its editor.
    ''')

    properties: SectionProperties = Field(None, description='''
        The value to customize the quantities and sub sections of the annotation section.
        The supported keys:
        `visible`: To determine the visible quantities and sub sections by their names<br/>
        `editable`: To render things visible but not editable, e.g. in inheritance situations<br/>
        `order`: # To order things, properties listed in that order first, then the rest<br/>
    ''')

    class Config:
        validate_assignment = True

    @validator('m_definition')
    def validate_component(cls, definition, values):  # pylint: disable=no-self-argument
        if not definition:
            return definition

        def assert_component(component, quantity_name, quantity_type, accepted_components):
            assert component in accepted_components, (
                f'The component {component} is not compatible with the quantity '
                f'{quantity_name} of the type {quantity_type}. '
                f'Accepted components: {", ".join(accepted_components)}.')

        component = values.get('component')
        if not component:
            return definition

        assert isinstance(definition, Quantity), 'Only quantities can be eln annotated with a component.'
        quantity = definition
        type_ = quantity.type
        name = quantity.name

        assert len(quantity.shape) <= 1, 'Only scalars or lists can be edited.'

        if isinstance(type_, type):
            if type_.__name__ == 'str':
                assert_component(component, name, type_.__name__, valid_eln_components['str'])
            elif type_.__name__ == 'bool':
                assert_component(component, name, type_.__name__, valid_eln_components['bool'])
            elif type_ in MTypes.num_python:
                assert_component(component, name, type_.__name__, valid_eln_components['number'])
            elif type_ in MTypes.num_numpy:
                assert_component(component, name, f'np.{type_.__name__}', valid_eln_components['number'])
            elif type_.__name__ == 'User':
                assert_component(component, name, type_.__name__, valid_eln_components['user'])
            elif type_.__name__ == 'Author':
                assert_component(component, name, type_.__name__, valid_eln_components['author'])

        elif type_ == Datetime:
            assert_component(component, name, type(type_).__name__, valid_eln_components['datetime'])

        elif isinstance(type_, MEnum):
            assert_component(component, name, type(type_).__name__, valid_eln_components['enum'])

        elif isinstance(type_, Reference):
            target_class = type_.target_section_def.section_cls
            if target_class.__name__ == 'User':
                assert_component(component, name, target_class.__name__, valid_eln_components['user'])
            elif target_class.__name__ == 'Author':
                assert_component(component, name, target_class.__name__, valid_eln_components['author'])
            else:
                assert_component(component, name, type(type_).__name__, valid_eln_components['reference'])

        return definition


class BrowserAdaptors(str, Enum):
    RawFileAdaptor = 'RawFileAdaptor'


class BrowserRenderValues(str, Enum):
    JsonValue = 'JsonValue'
    HtmlValue = 'HtmlValue'


class BrowserAnnotation(AnnotationModel):
    '''
    The `browser` annotation allows to specify if the processed data browser needs to
    display a quantity differently. It can be applied to quantities. For example

    ```python
        class Experiment(EntryData):
            description = Quantity(type=str, a_browser=dict(render_value='HtmlValue'))
    ```

    or in yaml

    ```yaml
    Experiment:
      quantities:
        description:
          type: str
          m_annotations:
            browser:
              render_value: HtmlValue
    ```
    '''

    adaptor: BrowserAdaptors = Field(None, description='''
      Allows to change the *Adaptor* implementation that is used to render the
      lane for this quantity. Possible values are:

      `RawFileAdaptor`: An adopter that is used to show files, including all file
      actions, like file preview.
    ''')
    render_value: BrowserRenderValues = Field(None, description='''
      Allows to change the *Component* used to render the value of the quantity.
      Possible values are:

      `HtmlValue`: Renders a string as HTML.<br/>
      `JsonValue`: Renders a dict or list in a collapsable tree.
    ''')


class TabularMode(str, Enum):
    row = 'row'
    column = 'column'
    root = 'root'
    entry = 'entry'


class TabularParserAnnotation(AnnotationModel):
    '''
    Instructs NOMAD to treat a string valued scalar quantity as a file path and
    interprets the contents of this file as tabular data. Supports both
    `.csv` and Excel files.
    '''

    comment: str = Field(None, description='''
        The character denoting the commented lines in `.csv` files. This is passed to
        pandas to parse the file. Has to be used to annotate the quantity that
        holds the path to the `.csv` or excel file.
    ''')
    sep: str = Field(None, description='''
        The character used to separate cells in a `.csv` file. This is passed to
        pandas to parse the file. Has to be used to annotate the quantity that
        holds the path to the `.csv` or excel file.
    ''')
    skiprows: int = Field(None, description='''
        Number of `.csv` file rows that are skipped. This is passed to
        pandas to parse the file. Has to be used to annotate the quantity that
        holds the path to the `.csv` or excel file.
    ''')
    separator: str = Field(None, description='An alias for `sep`')
    mode: TabularMode = Field(TabularMode.column, description='''
        Either `column`, `row`, `entry` or `root`. With `column` the whole column is mapped into a quantity
        (needs to be a list).
        With `row` each row (and its cells) are mapped into instances of a repeating
        sub section, where each section represents a row (quantities need to be scalars).
        With `entry` new entry is created and populated from each row (and its cells) where
        all quantities should remain to be scalars.
        Has to be used to annotate the quantity that
        holds the path to the `.csv` or excel file.
    ''')
    target_sub_section: List[str] = Field([], description='''
        A lists of paths to sub-sections of the annotation quantity's section. Each path is a
        `/` separated list of nested sub-sections. The targeted sub-sections, will be
        considered when mapping table columns to quantities.
        Has to be used to annotate the quantity that
        holds the path to the `.csv` or excel file.
    ''')


class TabularAnnotation(AnnotationModel):
    '''
    Allows to map a quantity to a row of a tabular data-file. Should only be used
    in conjunction with `tabular_parser`.
    '''

    name: str = Field(None, description='''
        The column name that should be mapped to the annotation quantity. Has to be
        the same string that is used in the header, i.e. first `.csv` line or first excel file `row`.
        For excel files with multiple sheets, the name can have the form `<sheet name>/<column name>`.
        Otherwise, only the first sheets is used. Has to be applied to the
        quantity that a column should be mapped to.
    ''')
    unit: str = Field(None, description='''
        The unit of the value in the file. Has to be compatible with the annotated quantity's
        unit. Will be used to automatically convert the value. If this is not defined,
        the values will not be converted. Has to be applied to the
        quantity that a column should be mapped to.
    ''')


class PlotAnnotation(AnnotationModel):
    '''
    This annotation can be used to add a plot to a section or quantity. Example:

    ```python
    class Evaporation(MSection):
        m_def = Section(a_plot={
            'label': 'Temperature and Pressure',
            'x': 'process_time',
            'y': ['./substrate_temperature', './chamber_pressure'],
            'config': {
                'editable': True,
                'scrollZoom': False
            }
        })
        time = Quantity(type=float, shape=['*'], unit='s')
        substrate_temperature = Quantity(type=float, shape=['*'], unit='K')
        chamber_pressure = Quantity(type=float, shape=['*'], unit='Pa')
    ```

    You can create multi-line plots by using lists of the properties `y` (and `x`).
    You either have multiple sets of `y`-values over a single set of `x`-values. Or
    you have pairs of `x` and `y` values. For this purpose the annotation properties
    `x` and `y` can reference a single quantity or a list of quantities.
    For repeating sub sections, the section instance can be selected with an index, e.g.
    "sub_section_name/2/parameter_name" or with a slice notation `start:stop` where
    negative values index from the end of the array, e.g.
    "sub_section_name/1:-5/parameter_name".

    The interactive examples of the plot annotations can be found
    [here](https://nomad-lab.eu/prod/v1/staging/gui/dev/plot).
    '''

    def __init__(self, *args, **kwargs):
        # pydantic does not seem to support multiple aliases per field
        super(PlotAnnotation, self).__init__(
            *args,
            x=kwargs.pop('x', None) or kwargs.pop('xAxis', None) or kwargs.pop('x_axis', None),
            y=kwargs.pop('y', None) or kwargs.pop('yAxis', None) or kwargs.pop('y_axis', None),
            **kwargs
        )

    label: str = Field(None, description='Is passed to plotly to define the label of the plot.')
    x: Union[List[str], str] = Field(..., description='''
        A path or list of paths to the x-axes values. Each path is a `/` separated
        list of sub-section and quantity names that leads from the annotation section
        to the quantity. Repeating sub sections are indexed between two `/`s with an
        integer or a slice `start:stop`.
    ''')
    y: Union[List[str], str] = Field(..., description='''
        A path or list of paths to the y-axes values. list of sub-section and quantity
        names that leads from the annotation section to the quantity. Repeating sub
        sections are indexed between two `/`s with an integer or a slice `start:stop`.
    ''')
    lines: List[dict] = Field(None, description='''
        A list of dicts passed as `traces` to plotly to configure the lines of the plot.
        See [https://plotly.com/javascript/reference/scatter/](https://plotly.com/javascript/reference/scatter/) for details.
    ''')
    layout: dict = Field(None, description='''
        A dict passed as `layout` to plotly to configure the plot layout.
        See [https://plotly.com/javascript/reference/layout/](https://plotly.com/javascript/reference/layout/) for details.
    ''')
    config: dict = Field(None, description='''
        A dict passed as `config` to plotly to configure the plot functionallity.
        See [https://plotly.com/javascript/configuration-options/](https://plotly.com/javascript/configuration-options/) for details.
    ''')

    @validator('y')
    def validate_y(cls, y, values):  # pylint: disable=no-self-argument
        x = values.get('x', [])
        if not isinstance(x, list):
            x = [x]

        if isinstance(x, list):
            assert len(x) == 1 or len(x) == len(y), strip(f'''
                You must use on set of x-values, or the amount x-quantities ({len(x)})
                has to match the amount of y-quantities ({len(y)}).
            ''')

        return y

    @validator('x', 'y')
    def validate_quantity_references(cls, value):  # pylint: disable=no-self-argument
        values = value if isinstance(value, list) else [value]
        for item in values:
            assert re.match(r'^(\.\/)?(\w+\/)*((\w+\/\-?\d*:\-?\d*)\/(\w+\/)*)*\w+$', item), f'{item} is not a valid quantity reference.'

        return value


AnnotationModel.m_registry['eln'] = ELNAnnotation
AnnotationModel.m_registry['browser'] = BrowserAnnotation
AnnotationModel.m_registry['tabular_parser'] = TabularParserAnnotation
AnnotationModel.m_registry['tabular'] = TabularAnnotation
AnnotationModel.m_registry['plot'] = PlotAnnotation
