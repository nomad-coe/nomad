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

from enum import Enum
import inspect
from typing import List, Dict, Union, Optional
from typing_extensions import Literal, Annotated
from pydantic import Field, root_validator

from .common import (
    ConfigBaseModel,
    Options,
    OptionsSingle,
    OptionsMulti,
    OptionsGlob,
    OptionsBase,
)


class UnitSystemUnit(ConfigBaseModel):
    definition: str = Field(
        description="""
        The unit definition. Can be a mathematical expression that combines
        several units, e.g. `(kg * m) / s^2`. You should only use units that are
        registered in the NOMAD unit registry (`nomad.units.ureg`).
    """
    )
    locked: Optional[bool] = Field(
        False,
        description='Whether the unit is locked in the unit system it is defined in.',
    )


dimensions = [
    # Base units
    'dimensionless',
    'length',
    'mass',
    'time',
    'current',
    'temperature',
    'luminosity',
    'luminous_flux',
    'substance',
    'angle',
    'information',
    # Derived units with specific name
    'force',
    'energy',
    'power',
    'pressure',
    'charge',
    'resistance',
    'conductance',
    'inductance',
    'magnetic_flux',
    'magnetic_field',
    'frequency',
    'luminance',
    'illuminance',
    'electric_potential',
    'capacitance',
    'activity',
]
dimension_list = '\n'.join([' - ' + str(dim) for dim in dimensions])


class UnitSystem(ConfigBaseModel):
    label: str = Field(
        description='Short, descriptive label used for this unit system.'
    )
    units: Optional[Dict[str, UnitSystemUnit]] = Field(
        description=f"""
        Contains a mapping from each dimension to a unit. If a unit is not
        specified for a dimension, the SI equivalent will be used by default.
        The following dimensions are available:
        {dimension_list}
    """
    )

    @root_validator(pre=True)
    def __validate(cls, values):  # pylint: disable=no-self-argument
        """Adds SI defaults for dimensions that are missing a unit."""
        units = values.get('units', {})
        from nomad.units import ureg
        from pint import UndefinedUnitError

        # Check that only supported dimensions and units are used
        for key in units.keys():
            if key not in dimensions:
                raise AssertionError(
                    f'Unsupported dimension "{key}" used in a unit system. The supported dimensions are: {dimensions}.'
                )

        # Fill missing units with SI defaults
        SI = {
            'dimensionless': 'dimensionless',
            'length': 'm',
            'mass': 'kg',
            'time': 's',
            'current': 'A',
            'temperature': 'K',
            'luminosity': 'cd',
            'luminous_flux': 'lm',
            'substance': 'mol',
            'angle': 'rad',
            'information': 'bit',
            'force': 'N',
            'energy': 'J',
            'power': 'W',
            'pressure': 'Pa',
            'charge': 'C',
            'resistance': 'Ω',
            'conductance': 'S',
            'inductance': 'H',
            'magnetic_flux': 'Wb',
            'magnetic_field': 'T',
            'frequency': 'Hz',
            'luminance': 'nit',
            'illuminance': 'lx',
            'electric_potential': 'V',
            'capacitance': 'F',
            'activity': 'kat',
        }
        for dimension in dimensions:
            if dimension not in units:
                units[dimension] = {'definition': SI[dimension]}

        # Check that units are available in registry, and thus also in the GUI.
        for value in units.values():
            definition = value['definition']
            try:
                ureg.Unit(definition)
            except UndefinedUnitError as e:
                raise AssertionError(
                    f'Unsupported unit "{definition}" used in a unit registry.'
                )

        values['units'] = units

        return values


class UnitSystems(OptionsSingle):
    """Controls the available unit systems."""

    options: Optional[Dict[str, UnitSystem]] = Field(
        description='Contains the available unit systems.'
    )


class Theme(ConfigBaseModel):
    """Theme and identity settings."""

    title: str = Field(description='Site name in the browser tab.')


class NORTHUI(ConfigBaseModel):
    """NORTH (NOMAD Remote Tools Hub) UI configuration."""

    enabled: bool = Field(
        True,
        description="""
        Whether the NORTH tools are available in the UI.
        The default value is read from the root-level NORTH configuration.
    """,
    )


class Card(ConfigBaseModel):
    """Definition for a card shown in the entry overview page."""

    error: str = Field(
        description='The error message to show if an error is encountered within the card.'
    )


class Cards(Options):
    """Contains the overview page card definitions and controls their visibility."""

    options: Optional[Dict[str, Card]] = Field(
        description='Contains the available card options.'
    )


class Entry(ConfigBaseModel):
    """Controls the entry visualization."""

    cards: Cards = Field(
        description='Controls the cards that are displayed on the entry overview page.'
    )


class Pagination(ConfigBaseModel):
    order_by: str = Field('upload_create_time', description='Field used for sorting.')
    order: str = Field('desc', description='Sorting order.')
    page_size: int = Field(20, description='Number of results on each page.')


class ModeEnum(str, Enum):
    STANDARD = 'standard'
    SCIENTIFIC = 'scientific'
    SEPARATORS = 'separators'


class Format(ConfigBaseModel):
    """Value formatting options."""

    decimals: int = Field(3, description='Number of decimals to show for numbers.')
    mode: ModeEnum = Field('standard', description='Display mode for numbers.')


class AlignEnum(str, Enum):
    LEFT = 'left'
    RIGHT = 'right'
    CENTER = 'center'


class Column(ConfigBaseModel):
    """Option for a column show in the search results."""

    label: Optional[str] = Field(
        description='Label shown in the header. Defaults to the quantity name.'
    )
    align: AlignEnum = Field(AlignEnum.LEFT, description='Alignment in the table.')
    unit: Optional[str] = Field(
        description="""
        Unit to convert to when displaying. If not given will be displayed in
        using the default unit in the active unit system.
    """
    )
    format: Optional[Format] = Field(
        description='Controls the formatting of the values.'
    )


class Columns(OptionsMulti):
    """
    Contains column definitions, controls their availability and specifies the default
    selection.
    """

    options: Optional[Dict[str, Column]] = Field(
        description="""
        All available column options. Note here that the key must correspond to a
        quantity path that exists in the metadata.
    """
    )


class RowActions(ConfigBaseModel):
    """Controls the visualization of row actions that are shown at the end of each row."""

    enabled: bool = Field(True, description='Whether to enable row actions. ')


class RowDetails(ConfigBaseModel):
    """
    Controls the visualization of row details that are shown upon pressing the row and
    contain basic details about the entry.
    """

    enabled: bool = Field(True, description='Whether to show row details.')


class RowSelection(ConfigBaseModel):
    """
    Controls the selection of rows. If enabled, rows can be selected and additional
    actions performed on them.
    """

    enabled: bool = Field(True, description='Whether to show the row selection.')


class Rows(ConfigBaseModel):
    """Controls the visualization of rows in the search results."""

    actions: RowActions
    details: RowDetails
    selection: RowSelection


class FilterMenuActionEnum(str, Enum):
    CHECKBOX = 'checkbox'


class FilterMenuAction(ConfigBaseModel):
    """Contains definition for an action in the filter menu."""

    type: FilterMenuActionEnum = Field(description='Action type.')
    label: str = Field(description='Label to show.')


class FilterMenuActionCheckbox(FilterMenuAction):
    """Contains definition for checkbox action in the filter menu."""

    quantity: str = Field(description='Targeted quantity')


class FilterMenuActions(Options):
    """Contains filter menu action definitions and controls their availability."""

    options: Optional[Dict[str, FilterMenuActionCheckbox]] = Field(
        description='Contains options for filter menu actions.'
    )


class FilterMenuSizeEnum(str, Enum):
    S = 's'
    M = 'm'
    L = 'l'
    XL = 'xl'


class FilterMenu(ConfigBaseModel):
    """Defines the layout and functionality for a filter menu."""

    label: Optional[str] = Field(description='Menu label to show in the UI.')
    level: Optional[int] = Field(0, description='Indentation level of the menu.')
    size: Optional[FilterMenuSizeEnum] = Field(
        FilterMenuSizeEnum.S, description='Width of the menu.'
    )
    actions: Optional[FilterMenuActions]


class FilterMenus(Options):
    """Contains filter menu definitions and controls their availability."""

    options: Optional[Dict[str, FilterMenu]] = Field(
        description='Contains the available filter menu options.'
    )


class Filters(OptionsGlob):
    """Controls the availability of filters in the app. Filters are pieces of
    (meta)info than can be queried in the search interface of the app, but also
    targeted in the rest of the app configuration. The `include` and `exlude`
    attributes can use glob syntax to target metainfo, e.g. `results.*` or
    `*.#myschema.schema.MySchema`.
    """

    include: Optional[List[str]] = Field(
        description="""
        List of included options. Supports glob/wildcard syntax.
    """
    )
    exclude: Optional[List[str]] = Field(
        description="""
        List of excluded options. Supports glob/wildcard syntax. Has higher precedence than include.
    """
    )


class SearchSyntaxes(ConfigBaseModel):
    """Controls the availability of different search syntaxes. These syntaxes
    determine how raw user input in e.g. the search bar is parsed into queries
    supported by the API.

    Currently you can only exclude items. By default, the following options are
    included:

     - `existence`: Used to query for the existence of a specific metainfo field in the data.
     - `equality`: Used to query for a specific value with exact match.
     - `range_bounded`: Queries values that are between two numerical limits, inclusive or exclusive.
     - `range_half_bounded`: Queries values that are above/below a numerical limit, inclusive or exclusive.
     - `free_text`: For inexact, free-text queries. Requires that a set of keywords has been filled in the entry.
    """

    exclude: Optional[List[str]] = Field(
        description="""
        List of excluded options.
    """
    )


class Layout(ConfigBaseModel):
    """Defines widget size and grid positioning for different breakpoints."""

    h: int = Field(description='Height in grid units')
    w: int = Field(description='Width in grid units.')
    x: int = Field(description='Horizontal start location in the grid.')
    y: int = Field(description='Vertical start location in the grid.')
    minH: Optional[int] = Field(3, description='Minimum height in grid units.')
    minW: Optional[int] = Field(3, description='Minimum width in grid units.')


class ScaleEnum(str, Enum):
    POW1 = 'linear'
    POW2 = '1/2'
    POW4 = '1/4'
    POW8 = '1/8'


class BreakpointEnum(str, Enum):
    SM = 'sm'
    MD = 'md'
    LG = 'lg'
    XL = 'xl'
    XXL = 'xxl'


class Axis(ConfigBaseModel):
    """Configuration for a plot axis."""

    title: Optional[str] = Field(description="""Custom title to show for the axis.""")
    unit: Optional[str] = Field(
        description="""Custom unit used for displaying the values."""
    )
    quantity: str = Field(
        description="""
        Path of the targeted quantity. Note that you can most of the features
        JMESPath syntax here to further specify a selection of values. This
        becomes especially useful when dealing with repeated sections or
        statistical values.
        """
    )


class Markers(ConfigBaseModel):
    """Configuration for plot markers."""

    color: Optional[Axis] = Field(
        description='Configures the information source and display options for the marker colors.'
    )


class Widget(ConfigBaseModel):
    """Common configuration for all widgets."""

    type: str = Field(description='Used to identify the widget type.')
    layout: Dict[BreakpointEnum, Layout] = Field(
        description="""
        Defines widget size and grid positioning for different breakpoints. The
        following breakpoints are supported: `sm`, `md`, `lg`, `xl` and `xxl`.
    """
    )


class WidgetTerms(Widget):
    """Terms widget configuration."""

    type: Literal['terms'] = Field(
        description='Set as `terms` to get this widget type.'
    )
    quantity: str = Field(description='Targeted quantity.')
    scale: ScaleEnum = Field(description='Statistics scaling.')
    showinput: bool = Field(True, description='Whether to show text input field.')


class WidgetHistogram(Widget):
    """Histogram widget configuration."""

    type: Literal['histogram'] = Field(
        description='Set as `histogram` to get this widget type.'
    )
    quantity: str = Field(description='Targeted quantity.')
    scale: ScaleEnum = Field(description='Statistics scaling.')
    autorange: bool = Field(
        True,
        description='Whether to automatically set the range according to the data limits.',
    )
    showinput: bool = Field(
        True,
        description='Whether to show input text fields for minimum and maximum value.',
    )
    nbins: int = Field(
        description="""
        Maximum number of histogram bins. Notice that the actual number of bins
        may be smaller if there are fewer data items available.
        """
    )


class WidgetPeriodicTable(Widget):
    """Periodic table widget configuration."""

    type: Literal['periodictable'] = Field(
        description='Set as `periodictable` to get this widget type.'
    )
    quantity: str = Field(description='Targeted quantity.')
    scale: ScaleEnum = Field(description='Statistics scaling.')


class WidgetScatterPlot(Widget):
    """Scatter plot widget configuration."""

    type: Literal['scatterplot'] = Field(
        description='Set as `scatterplot` to get this widget type.'
    )
    x: Union[Axis, str] = Field(
        description='Configures the information source and display options for the x-axis.'
    )
    y: Union[Axis, str] = Field(
        description='Configures the information source and display options for the y-axis.'
    )
    markers: Optional[Markers] = Field(
        description='Configures the information source and display options for the markers.'
    )
    color: Optional[str] = Field(
        description="""
        Quantity used for coloring points. Note that this field is deprecated
        and `markers` should be used instead.
        """
    )
    size: int = Field(
        1000,
        description="""
        Maximum number of entries to fetch. Notice that the actual number may be
        more of less, depending on how many entries exist and how many of the
        requested values each entry contains.
        """,
    )
    autorange: bool = Field(
        True,
        description='Whether to automatically set the range according to the data limits.',
    )

    @root_validator(pre=True)
    def __validate(cls, values):
        """Ensures backwards compatibility of x, y, and color."""
        color = values.get('color')
        if color is not None:
            values['markers'] = {'color': {'quantity': color}}
            del values['color']
        x = values.get('x')
        if isinstance(x, str):
            values['x'] = {'quantity': x}
        y = values.get('y')
        if isinstance(y, str):
            values['y'] = {'quantity': y}
        return values


# The 'discriminated union' feature of Pydantic is used here:
# https://docs.pydantic.dev/usage/types/#discriminated-unions-aka-tagged-unions
WidgetAnnotated = Annotated[
    Union[WidgetTerms, WidgetHistogram, WidgetScatterPlot, WidgetPeriodicTable],
    Field(discriminator='type'),
]


class Dashboard(ConfigBaseModel):
    """Dashboard configuration."""

    widgets: List[WidgetAnnotated] = Field(
        description='List of widgets contained in the dashboard.'
    )  # type: ignore


class ResourceEnum(str, Enum):
    ENTRIES = 'entries'
    MATERIALS = 'materials'


class App(ConfigBaseModel):
    """Defines the layout and functionality for an App."""

    label: str = Field(description='Name of the App.')
    path: str = Field(description='Path used in the browser address bar.')
    resource: ResourceEnum = Field('entries', description='Targeted resource.')
    breadcrumb: Optional[str] = Field(
        description='Name displayed in the breadcrumb, by default the label will be used.'
    )
    category: str = Field(
        description='Category used to organize Apps in the explore menu.'
    )
    description: Optional[str] = Field(description='Short description of the App.')
    readme: Optional[str] = Field(
        description='Longer description of the App that can also use markdown.'
    )
    pagination: Pagination = Field(
        Pagination(), description='Default result pagination.'
    )
    columns: Columns = Field(
        description='Controls the columns shown in the results table.'
    )
    rows: Optional[Rows] = Field(
        Rows(
            actions=RowActions(enabled=True),
            details=RowDetails(enabled=True),
            selection=RowSelection(enabled=True),
        ),
        description='Controls the display of entry rows in the results table.',
    )
    filter_menus: FilterMenus = Field(
        description='Filter menus displayed on the left side of the screen.'
    )
    filters: Optional[Filters] = Field(
        Filters(exclude=['mainfile', 'entry_name', 'combine']),
        description='Controls the filters that are available in this app.',
    )
    dashboard: Optional[Dashboard] = Field(description='Default dashboard layout.')
    filters_locked: Optional[dict] = Field(
        description="""
        Fixed query object that is applied for this search context. This filter
        will always be active for this context and will not be displayed to the
        user by default.
        """
    )
    search_syntaxes: Optional[SearchSyntaxes] = Field(
        description='Controls which types of search syntax are available.'
    )


class Apps(Options):
    """Contains App definitions and controls their availability."""

    options: Optional[Dict[str, App]] = Field(
        description='Contains the available app options.'
    )


class ExampleUploads(OptionsBase):
    """Controls the availability of example uploads."""


class UI(ConfigBaseModel):
    """Used to customize the user interface."""

    app_base: str = Field(None, description='This is automatically set.')
    north_base: str = Field(None, description='This is automatically set.')
    theme: Theme = Field(None, description='Controls the site theme and identity.')
    unit_systems: UnitSystems = Field(
        None, description='Controls the available unit systems.'
    )
    entry: Entry = Field(None, description='Controls the entry visualization.')
    apps: Apps = Field(None, description='Contains the App definitions.')
    north: NORTHUI = Field(
        NORTHUI(), description='NORTH (NOMAD Remote Tools Hub) UI configuration.'
    )
    example_uploads: ExampleUploads = Field(
        ExampleUploads(), description='Controls the available example uploads.'
    )
