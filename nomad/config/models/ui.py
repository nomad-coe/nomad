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


class ScaleEnum(str, Enum):
    LINEAR = 'linear'
    LOG = 'log'
    # TODO: The following should possibly be deprecated.
    POW1 = 'linear'
    POW2 = '1/2'
    POW4 = '1/4'
    POW8 = '1/8'


class ScaleEnumPlot(str, Enum):
    LINEAR = 'linear'
    LOG = 'log'


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
    def _validate(cls, values):  # pylint: disable=no-self-argument
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
    DATE = 'date'
    TIME = 'time'


class Format(ConfigBaseModel):
    """Value formatting options."""

    decimals: int = Field(3, description='Number of decimals to show for numbers.')
    mode: ModeEnum = Field(ModeEnum.SCIENTIFIC, description='Display mode for numbers.')


class AlignEnum(str, Enum):
    LEFT = 'left'
    RIGHT = 'right'
    CENTER = 'center'


class Column(ConfigBaseModel):
    """Column show in the search results table. With `quantity` you may target a
    specific part of the data to be shown. Note that the use of JMESPath is
    supported here, and you can e.g. do the following:

     - Show first value from a repeating subsection: `repeating_section[0].quantity`
     - Show slice of values from a repeating subsection: `repeating_section[1:2].quantity`
     - Show all values from a repeating subsection: `repeating_section[*].quantity`
     - Show minimum value from a repeating section: `min(repeating_section[*].quantity)`
     - Show instance that matches a criterion: `repeating_section[?label=='target'].quantity`
    """

    search_quantity: Optional[str] = Field(
        description="""
        Path of the targeted quantity. Note that you can most of the features
        JMESPath syntax here to further specify a selection of values. This
        becomes especially useful when dealing with repeated sections or
        statistical values.
        """
    )
    quantity: Optional[str] = Field(
        deprecated='The "quantity" field is deprecated, use "search_quantity" instead.'
    )
    selected: bool = Field(
        False, description="""Is this column initially selected to be shown."""
    )
    title: Optional[str] = Field(
        description='Label shown in the header. Defaults to the quantity name.'
    )
    label: Optional[str] = Field(description='Alias for title.')
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

    @root_validator(pre=True)
    def _validate(cls, values):
        # Backwards compatibility for quantity.
        quantity = values.get('quantity')
        search_quantity = values.get('search_quantity')
        if quantity is not None and search_quantity is None:
            values['search_quantity'] = quantity
            del values['quantity']

        # Backwards compatibility for label
        label = values.get('label')
        title = values.get('title')
        if label and not title:
            values['title'] = label
        return values


class Columns(OptionsMulti):
    """
    Contains column definitions, controls their availability and specifies the default
    selection.
    """

    options: Optional[Dict[str, Column]] = Field(
        description="""
        All available column options. Note here that the key must correspond to a
        quantity path that exists in the metadata. The key
    """
    )


class RowAction(ConfigBaseModel):
    """Common configuration for all row actions."""

    description: Optional[str] = Field(
        description="""Description of the action shown to the user."""
    )
    type: str = Field(description='Used to identify the action type.')


class RowActionURL(RowAction):
    """Action that will open an external link read from the archive."""

    path: str = Field(
        description="""JMESPath pointing to a path in the archive that contains the URL."""
    )
    type: Literal['url'] = Field(
        'url', description='Set as `url` to get this widget type.'
    )


class RowActions(Options):
    """Controls the visualization of row actions that are shown at the end of each row."""

    enabled: bool = Field(True, description='Whether to enable row actions.')
    options: Optional[Dict[str, RowActionURL]] = Field(
        description="""
        All available row actions.
    """
    )


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


# Deprecated
class FilterMenuActionEnum(str, Enum):
    CHECKBOX = 'checkbox'


# Deprecated
class FilterMenuAction(ConfigBaseModel):
    """Contains definition for an action in the filter menu."""

    type: FilterMenuActionEnum = Field(description='Action type.')
    label: str = Field(description='Label to show.')


# Deprecated
class FilterMenuActionCheckbox(FilterMenuAction):
    """Contains definition for checkbox action in the filter menu."""

    quantity: str = Field(description='Targeted quantity')


# Deprecated
class FilterMenuActions(Options):
    """Contains filter menu action definitions and controls their availability."""

    options: Optional[Dict[str, FilterMenuActionCheckbox]] = Field(
        description='Contains options for filter menu actions.'
    )


# Deprecated
class FilterMenuSizeEnum(str, Enum):
    S = 's'
    M = 'm'
    L = 'l'
    XL = 'xl'


# Deprecated
class FilterMenu(ConfigBaseModel):
    """Defines the layout and functionality for a filter menu."""

    label: Optional[str] = Field(description='Menu label to show in the UI.')
    level: Optional[int] = Field(0, description='Indentation level of the menu.')
    size: Optional[FilterMenuSizeEnum] = Field(
        FilterMenuSizeEnum.S, description='Width of the menu.'
    )
    actions: Optional[FilterMenuActions]


# Deprecated
class FilterMenus(Options):
    """Contains filter menu definitions and controls their availability."""

    options: Optional[Dict[str, FilterMenu]] = Field(
        description='Contains the available filter menu options.'
    )


# NOTE: Once the old power scaling options (1/2, 1/4, 1/8) are deprecated, the
# axis models here can be simplified.
class AxisScale(ConfigBaseModel):
    """Basic configuration for a plot axis."""

    scale: Optional[ScaleEnum] = Field(
        ScaleEnum.LINEAR,
        description="""Defines the axis scaling. Defaults to linear scaling.""",
    )


class AxisQuantity(ConfigBaseModel):
    """Configuration for a plot axis."""

    title: Optional[str] = Field(description="""Custom title to show for the axis.""")
    unit: Optional[str] = Field(
        description="""Custom unit used for displaying the values."""
    )
    quantity: Optional[str] = Field(
        deprecated='The "quantity" field is deprecated, use "search_quantity" instead.'
    )
    search_quantity: str = Field(
        description="""
        Path of the targeted search quantity. Note that you can most of the features
        JMESPath syntax here to further specify a selection of values. This
        becomes especially useful when dealing with repeated sections or
        statistical values.
        """
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        # Backwards compatibility for quantity.
        quantity = values.get('quantity')
        search_quantity = values.get('search_quantity')
        if quantity is not None and search_quantity is None:
            values['search_quantity'] = quantity
            del values['quantity']

        return values


class Axis(AxisScale, AxisQuantity):
    """Configuration for a plot axis with limited scaling options."""


class TermsBase(ConfigBaseModel):
    """Base model for configuring terms components."""

    quantity: Optional[str] = Field(
        deprecated='The "quantity" field is deprecated, use "search_quantity" instead.'
    )
    search_quantity: str = Field(description='The targeted search quantity.')
    type: Literal['terms'] = Field(
        description='Set as `terms` to get this type.',
    )
    scale: ScaleEnum = Field(ScaleEnum.LINEAR, description='Statistics scaling.')
    show_input: bool = Field(True, description='Whether to show text input field.')
    showinput: Optional[bool] = Field(
        deprecated='The "showinput" field is deprecated, use "show_input" instead.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'terms'

        # Backwards compatibility for showinput.
        showinput = values.get('showinput')
        if showinput is not None:
            values['show_input'] = showinput

        # Backwards compatibility for quantity.
        quantity = values.get('quantity')
        search_quantity = values.get('search_quantity')
        if quantity is not None and search_quantity is None:
            values['search_quantity'] = quantity
            del values['quantity']

        return values


class HistogramBase(ConfigBaseModel):
    """Base model for configuring histogram components."""

    type: Literal['histogram'] = Field(
        description='Set as `histogram` to get this widget type.'
    )
    quantity: Optional[str] = Field(
        deprecated='The "quantity" field is deprecated, use "x.search_quantity" instead.'
    )
    scale: Optional[ScaleEnum] = Field(
        deprecated='The "scale" field is deprecated, use "y.scale" instead.'
    )
    show_input: bool = Field(True, description='Whether to show text input field.')
    showinput: Optional[bool] = Field(
        deprecated='The "showinput" field is deprecated, use "show_input" instead.'
    )

    x: Union[Axis, str] = Field(
        description='Configures the information source and display options for the x-axis.'
    )
    y: Union[AxisScale, str] = Field(
        description='Configures the information source and display options for the y-axis.'
    )
    autorange: bool = Field(
        False,
        description='Whether to automatically set the range according to the data limits.',
    )
    n_bins: Optional[int] = Field(
        description="""
        Maximum number of histogram bins. Notice that the actual number of bins
        may be smaller if there are fewer data items available.
        """
    )
    nbins: Optional[int] = Field(
        deprecated='The "nbins" field is deprecated, use "n_bins" instead.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'histogram'

        # Backwards compatibility for nbins."""
        nbins = values.get('nbins')
        if nbins is not None:
            values['n_bins'] = nbins

        # Backwards compatibility for showinput."""
        showinput = values.get('showinput')
        if showinput is not None:
            values['show_input'] = showinput

        # x backwards compatibility
        x = values.get('x', {})
        if isinstance(x, str):
            x = {'search_quantity': x}
        if isinstance(x, dict):
            quantity = values.get('quantity')
            if quantity and not x.get('search_quantity'):
                x['search_quantity'] = quantity
                del values['quantity']
            values['x'] = x

        # y backwards compatibility
        y = values.get('y', {})
        if isinstance(y, dict):
            scale = values.get('scale')
            if scale:
                y['scale'] = scale
                del values['scale']
            values['y'] = y

        return values


class PeriodicTableBase(ConfigBaseModel):
    """Base model for configuring periodic table components."""

    type: Literal['periodic_table'] = Field(
        description='Set as `periodic_table` to get this widget type.'
    )
    quantity: Optional[str] = Field(
        deprecated='The "quantity" field is deprecated, use "search_quantity" instead.'
    )
    search_quantity: str = Field(description='The targeted search quantity.')
    scale: Optional[ScaleEnum] = Field(
        ScaleEnum.LINEAR, description='Statistics scaling.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'periodic_table'

        # Backwards compatibility for quantity.
        quantity = values.get('quantity')
        search_quantity = values.get('search_quantity')
        if quantity is not None and search_quantity is None:
            values['search_quantity'] = quantity
            del values['quantity']

        return values


class MenuSizeEnum(str, Enum):
    XS = 'xs'
    SM = 'sm'
    MD = 'md'
    LG = 'lg'
    XL = 'xl'
    XXL = 'xxl'


class MenuItem(ConfigBaseModel):
    width: int = Field(
        12,
        description='Width of the item, 12 means maximum width. Note that the menu size can be changed.',
    )
    show_header: bool = Field(True, description='Whether to show the header.')
    title: Optional[str] = Field(description='Custom item title.')


class MenuItemOption(ConfigBaseModel):
    """Represents an option shown for a filter."""

    label: Optional[str] = Field(description='The label to show for this option.')
    description: Optional[str] = Field(
        description='Detailed description for this option.'
    )


class MenuItemTerms(MenuItem, TermsBase):
    """Menu item that shows a list of text values from e.g. `str` or `MEnum`
    quantities.
    """

    options: Optional[Union[int, Dict[str, MenuItemOption]]] = Field(
        description="""
        Used to control the displayed options:

         - If not specified, sensible default options are shown based on the
           definition. For enum fields all of the defined options are shown,
           whereas for generic string fields the top 5 options are shown.

         - If a number is specified, that many options are dynamically fetched
           in order of occurrence. Set to 0 to completely disable options.

         - If a dictionary of str + MenuItemOption pairs is given, only these
           options will be shown.
        """
    )
    n_columns: int = Field(
        1,
        description='The number of columns to use when displaying the options.',
    )
    sort_static: bool = Field(
        True,
        description="""
        Whether to sort static options by their occurrence in the data. Options
        are static if they are read from the enum options of the field or if
        they are explicitly given as a dictionary in 'options'.
        """,
    )
    show_statistics: bool = Field(
        True, description='Whether to show statistics for the options.'
    )


class MenuItemHistogram(MenuItem, HistogramBase):
    """Menu item that shows a histogram for numerical or timestamp quantities."""

    show_statistics: bool = Field(
        True, description='Whether to show the full histogram, or just a range slider.'
    )


class MenuItemPeriodicTable(MenuItem, PeriodicTableBase):
    """Menu item that shows a periodic table built from values stored into a
    text quantity.
    """

    show_statistics: bool = Field(
        True, description='Whether to show statistics for the options.'
    )


class MenuItemVisibility(MenuItem):
    """Menu item that shows a radio button that can be used to change the visiblity."""

    type: Literal['visibility'] = Field(
        description='Set as `visibility` to get this menu item type.',
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'visibility'
        return values


class MenuItemDefinitions(MenuItem):
    """Menu item that shows a tree for filtering data by the presence of definitions."""

    type: Literal['definitions'] = Field(
        description='Set as `definitions` to get this menu item type.',
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'definitions'
        return values


class MenuItemOptimade(MenuItem):
    """Menu item that shows a dialog for entering OPTIMADE queries."""

    type: Literal['optimade'] = Field(
        description='Set as `optimade` to get this menu item type.',
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'optimade'
        return values


class MenuItemCustomQuantities(MenuItem):
    """Menu item that shows a search dialog for filtering by custom quantities
    coming from all different custom schemas, including YAML and Python schemas.
    Will only show quantities that have been populated in the data.
    """

    type: Literal['custom_quantities'] = Field(
        description='Set as `custom_quantities` to get this menu item type.',
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'custom_quantities'
        return values


# The 'discriminated union' feature of Pydantic is used here:
# https://docs.pydantic.dev/usage/types/#discriminated-unions-aka-tagged-unions
MenuItemTypeNested = Annotated[
    Union[
        MenuItemTerms,
        MenuItemHistogram,
        MenuItemPeriodicTable,
        MenuItemVisibility,
        MenuItemDefinitions,
        MenuItemOptimade,
        MenuItemCustomQuantities,
    ],
    Field(discriminator='type'),
]


class MenuItemNestedObject(MenuItem):
    """Menu item that can be used to wrap several subitems into a nested object.
    By wrapping items with this class the query for them is performed as an
    Elasticsearch nested query:
    https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-nested-query.html.
    Note that you cannot yet use nested queries for search quantities
    originating from custom schemas.
    """

    type: Literal['nested_object'] = Field(
        description='Set as `nested_object` to get this menu item type.',
    )
    path: str = Field(
        description='Path of the nested object. Typically a section name.'
    )
    items: Optional[List[MenuItemTypeNested]] = Field(
        description='Items that are grouped by this nested object.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'nested_object'
        return values


# The 'discriminated union' feature of Pydantic is used here:
# https://docs.pydantic.dev/usage/types/#discriminated-unions-aka-tagged-unions
MenuItemType = Annotated[
    Union[
        MenuItemTerms,
        MenuItemHistogram,
        MenuItemPeriodicTable,
        MenuItemNestedObject,
        MenuItemVisibility,
        MenuItemDefinitions,
        MenuItemOptimade,
        MenuItemCustomQuantities,
        'Menu',
    ],
    Field(discriminator='type'),
]


class Menu(MenuItem):
    """Defines a menu that is shown on the left side of the search interface.
    Menus have a controllable width, and contains items. Items in the menu are
    displayed on a 12-based grid and you can control the width of each item by
    using the `width` field. You can also nest menus within each other.
    """

    type: Literal['menu'] = Field(
        description='Set as `nested_object` to get this menu item type.',
    )
    size: Optional[Union[MenuSizeEnum, str]] = Field(
        MenuSizeEnum.SM,
        description="""
        Size of the menu. Either use presets as defined by MenuSizeEnum,
        or then provide valid CSS widths.
        """,
    )
    indentation: Optional[int] = Field(0, description='Indentation level for the menu.')
    items: Optional[List[MenuItemType]] = Field(
        description='List of items in the menu.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'menu'
        return values


class SearchQuantities(OptionsGlob):
    """Controls the quantities that are available in the search interface.
    Search quantities correspond to pieces of information that can be queried in
    the search interface of the app, but also targeted in the rest of the app
    configuration. You can load quantities from custom schemas as search
    quantities, but note that not all quantities will be loaded: only scalar
    values are supported at the moment. The `include` and `exlude` attributes
    can use glob syntax to target metainfo, e.g. `results.*` or
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


class Filters(SearchQuantities):
    """Alias for SearchQuantities."""


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


class BreakpointEnum(str, Enum):
    SM = 'sm'
    MD = 'md'
    LG = 'lg'
    XL = 'xl'
    XXL = 'xxl'


class AxisLimitedScale(AxisQuantity):
    """Configuration for a plot axis with limited scaling options."""

    scale: Optional[ScaleEnumPlot] = Field(
        ScaleEnumPlot.LINEAR,
        description="""Defines the axis scaling. Defaults to linear scaling.""",
    )


class Markers(ConfigBaseModel):
    """Configuration for plot markers."""

    color: Optional[Axis] = Field(
        description='Configures the information source and display options for the marker colors.'
    )


class Widget(ConfigBaseModel):
    """Common configuration for all widgets."""

    title: Optional[str] = Field(
        description='Custom widget title. If not specified, a widget-specific default title is used.'
    )
    type: str = Field(description='Used to identify the widget type.')
    layout: Dict[BreakpointEnum, Layout] = Field(
        description="""
        Defines widget size and grid positioning for different breakpoints. The
        following breakpoints are supported: `sm`, `md`, `lg`, `xl` and `xxl`.
    """
    )


class WidgetTerms(Widget, TermsBase):
    """Terms widget configuration."""

    type: Literal['terms'] = Field(
        description='Set as `terms` to get this type.',
    )


class WidgetHistogram(Widget, HistogramBase):
    """Histogram widget configuration."""

    type: Literal['histogram'] = Field(
        description='Set as `histogram` to get this type.',
    )


class WidgetPeriodicTable(Widget, PeriodicTableBase):
    """Periodic table widget configuration."""

    type: Literal['periodic_table'] = Field(
        description='Set as `periodic_table` to get this type.',
    )


class WidgetPeriodicTableDeprecated(WidgetPeriodicTable):
    """Deprecated copy of WidgetPeriodicTable with a misspelled type."""

    type: Literal['periodictable'] = Field(  # type: ignore[assignment]
        description='Set as `periodictable` to get this widget type.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'periodictable'

        # Backwards compatibility for quantity.
        quantity = values.get('quantity')
        search_quantity = values.get('search_quantity')
        if quantity is not None and search_quantity is None:
            values['search_quantity'] = quantity
            del values['quantity']

        return values


class WidgetScatterPlot(Widget):
    """Scatter plot widget configuration."""

    type: Literal['scatter_plot'] = Field(
        description='Set as `scatter_plot` to get this widget type.'
    )
    x: Union[AxisLimitedScale, str] = Field(
        description='Configures the information source and display options for the x-axis.'
    )
    y: Union[AxisLimitedScale, str] = Field(
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
        more or less, depending on how many entries exist and how many of the
        requested values each entry contains.
        """,
    )
    autorange: bool = Field(
        True,
        description='Whether to automatically set the range according to the data limits.',
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'scatter_plot'

        # color backwards compatibility
        color = values.get('color')
        if color is not None:
            values['markers'] = {'color': {'search_quantity': color}}
            del values['color']

        # x backwards compatibility
        x = values.get('x')
        if isinstance(x, str):
            values['x'] = {'search_quantity': x}

        # y backwards compatibility
        y = values.get('y')
        if isinstance(y, str):
            values['y'] = {'search_quantity': y}
        return values


class WidgetScatterPlotDeprecated(WidgetScatterPlot):
    """Deprecated copy of WidgetScatterPlot with a misspelled type."""

    type: Literal['scatterplot'] = Field(  # type: ignore[assignment]
        description='Set as `scatterplot` to get this widget type.'
    )

    @root_validator(pre=True)
    def _validate(cls, values):
        values['type'] = 'scatterplot'
        return values


# The 'discriminated union' feature of Pydantic is used here:
# https://docs.pydantic.dev/usage/types/#discriminated-unions-aka-tagged-unions
WidgetAnnotated = Annotated[
    Union[
        WidgetTerms,
        WidgetHistogram,
        WidgetScatterPlot,
        WidgetScatterPlotDeprecated,
        WidgetPeriodicTable,
        WidgetPeriodicTableDeprecated,
    ],
    Field(discriminator='type'),
]


class Dashboard(ConfigBaseModel):
    """Dashboard configuration."""

    widgets: List[WidgetAnnotated] = Field(
        description='List of widgets contained in the dashboard.'
    )


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
    columns: Optional[List[Column]] = Field(
        description='List of columns for the results table.'
    )
    rows: Optional[Rows] = Field(
        Rows(
            actions=RowActions(enabled=True),
            details=RowDetails(enabled=True),
            selection=RowSelection(enabled=True),
        ),
        description='Controls the display of entry rows in the results table.',
    )
    menu: Optional[Menu] = Field(
        description='Filter menu displayed on the left side of the screen.'
    )
    filter_menus: Optional[FilterMenus] = Field(
        deprecated='The "filter_menus" field is deprecated, use "menu" instead.'
    )
    filters: Optional[Filters] = Field(
        deprecated='The "filters" field is deprecated, use "search_quantities" instead.'
    )
    search_quantities: Optional[SearchQuantities] = Field(
        SearchQuantities(exclude=['mainfile', 'entry_name', 'combine']),
        description='Controls the quantities that are available for search in this app.',
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

    @root_validator(pre=True)
    def _validate(cls, values):
        # Backwards compatibility for columns
        columns = values.get('columns')
        if isinstance(columns, Columns):
            columns = columns.dict()
        if isinstance(columns, dict):
            options = columns.get('options') or {}
            keys = list(options.keys())
            include = columns.get('include') or keys
            exclude = columns.get('exclude') or []
            selected = columns.get('selected') or []
            new_columns = []
            for key in include:
                if key not in exclude:
                    column = options[key]
                    column['quantity'] = key
                    if key in selected:
                        column['selected'] = True
                    new_columns.append(column)
            values['columns'] = new_columns

        # Backwards compatibility for Filters
        filters = values.get('filters')
        if filters and not values.get('search_quantities'):
            values['search_quantities'] = filters

        # Backwards compatibility for FilterMenus
        filter_menus = values.get('filter_menus')
        if isinstance(filter_menus, FilterMenus):
            filter_menus = filter_menus.dict()
        options = filter_menus.get('options') if filter_menus else None
        menus = values.get('menus')
        if options and not menus:
            items = []
            for key, value in options.items():
                menu = {
                    'material': Menu(),
                    'elements': Menu(
                        items=[
                            MenuItemPeriodicTable(
                                search_quantity='results.material.elements',
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.chemical_formula_hill',
                                width=6,
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.chemical_formula_iupac',
                                width=6,
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.chemical_formula_reduced',
                                width=6,
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.chemical_formula_anonymous',
                                width=6,
                                options=0,
                            ),
                            MenuItemHistogram(
                                x='results.material.n_elements',
                            ),
                        ]
                    ),
                    'structure': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.material.structural_type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.bravais_lattice',
                                n_columns=2,
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.crystal_system',
                                n_columns=2,
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.space_group_symbol',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.structure_name',
                                options=5,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.strukturbericht_designation',
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.point_group',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.hall_symbol',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.symmetry.prototype_aflow_id',
                                options=0,
                            ),
                        ]
                    ),
                    'method': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.simulation.program_name',
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.program_version',
                                options=0,
                            ),
                        ]
                    ),
                    'precision': Menu(
                        items=[
                            MenuItemHistogram(
                                x='results.method.simulation.precision.k_line_density',
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.precision.native_tier',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.precision.basis_set',
                                options=5,
                            ),
                            MenuItemHistogram(
                                x=Axis(
                                    search_quantity='results.method.simulation.precision.planewave_cutoff',
                                )
                            ),
                            MenuItemHistogram(
                                x='results.method.simulation.precision.apw_cutoff',
                            ),
                        ]
                    ),
                    'dft': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.method_name',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'DFT': MenuItemOption(label='Search DFT entries')
                                },
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dft.xc_functional_type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dft.xc_functional_names',
                            ),
                            MenuItemHistogram(
                                x='results.method.simulation.dft.exact_exchange_mixing_factor',
                            ),
                            MenuItemHistogram(
                                x='results.method.simulation.dft.hubbard_kanamori_model.u_effective',
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dft.core_electron_treatment',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dft.relativity_method',
                                show_input=False,
                            ),
                        ]
                    ),
                    'tb': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.method_name',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'TB': MenuItemOption(label='Search TB entries')
                                },
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.tb.type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.tb.localization_type',
                                show_input=False,
                            ),
                        ]
                    ),
                    'gw': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.method_name',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'GW': MenuItemOption(label='Search GW entries')
                                },
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.gw.type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.gw.starting_point_type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.gw.basis_set_type',
                                show_input=False,
                            ),
                        ]
                    ),
                    'bse': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.method_name',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'BSE': MenuItemOption(label='Search BSE entries')
                                },
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.bse.type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.bse.solver',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.bse.starting_point_type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.bse.starting_point_type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.bse.basis_set_type',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.bse.gw_type',
                                show_input=False,
                            ),
                        ]
                    ),
                    'dmft': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.method_name',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'DMFT': MenuItemOption(label='Search DMFT entries')
                                },
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dmft.impurity_solver_type',
                                n_columns=2,
                                show_input=False,
                            ),
                            MenuItemHistogram(
                                x='results.method.simulation.dmft.inverse_temperature',
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dmft.magnetic_state',
                                show_input=False,
                            ),
                            MenuItemHistogram(
                                x='results.method.simulation.dmft.u',
                            ),
                            MenuItemHistogram(
                                x='results.method.simulation.dmft.jh',
                            ),
                            MenuItemTerms(
                                search_quantity='results.method.simulation.dmft.analytical_continuation',
                                show_input=False,
                            ),
                        ]
                    ),
                    'eels': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.method.method_name',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'EELS': MenuItemOption(label='Search EELS entries')
                                },
                            ),
                            MenuItemNestedObject(
                                path='results.properties.spectroscopic.spectra.provenance.eels',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.spectroscopic.spectra.provenance.eels.detector_type'
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.spectroscopic.spectra.provenance.eels.resolution',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.spectroscopic.spectra.provenance.eels.min_energy',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.spectroscopic.spectra.provenance.eels.max_energy',
                                    ),
                                ],
                            ),
                        ]
                    ),
                    'workflow': Menu(),
                    'molecular_dynamics': Menu(
                        items=[
                            MenuItemNestedObject(
                                path='results.properties.thermodynamic.trajectory',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.thermodynamic.trajectory.available_properties',
                                        show_input=False,
                                        options=4,
                                    ),
                                    MenuItemTerms(
                                        search_quantity='results.properties.thermodynamic.trajectory.provenance.molecular_dynamics.ensemble_type',
                                        show_input=False,
                                        options=2,
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.thermodynamic.trajectory.provenance.molecular_dynamics.time_step',
                                    ),
                                ],
                            ),
                        ]
                    ),
                    'geometry_optimization': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.properties.available_properties',
                                show_header=False,
                                show_input=False,
                                show_statistics=False,
                                options={
                                    'geometry_optimization': MenuItemOption(
                                        label='Search geometry optimization entries'
                                    )
                                },
                            ),
                            MenuItemNestedObject(
                                path='results.properties.geometry_optimization',
                                items=[
                                    MenuItemHistogram(
                                        x='results.properties.geometry_optimization.final_energy_difference',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.geometry_optimization.final_force_maximum',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.geometry_optimization.final_displacement_maximum',
                                    ),
                                ],
                            ),
                        ]
                    ),
                    'properties': Menu(),
                    'electronic': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='electronic_properties',
                                show_input=False,
                            ),
                            MenuItemNestedObject(
                                path='results.properties.electronic.band_structure_electronic.band_gap',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.electronic.band_structure_electronic.band_gap.type',
                                        options=2,
                                        show_input=False,
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.electronic.band_structure_electronic.band_gap.value',
                                    ),
                                ],
                            ),
                            MenuItemNestedObject(
                                path='results.properties.electronic.band_structure_electronic',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.electronic.band_structure_electronic.spin_polarized',
                                        show_input=False,
                                    ),
                                ],
                            ),
                            MenuItemNestedObject(
                                path='results.properties.electronic.dos_electronic',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.electronic.dos_electronic.spin_polarized',
                                        show_input=False,
                                    ),
                                ],
                            ),
                        ]
                    ),
                    'vibrational': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='vibrational_properties',
                                show_input=False,
                            ),
                        ]
                    ),
                    'mechanical': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='mechanical_properties',
                                show_input=False,
                            ),
                            MenuItemNestedObject(
                                path='results.properties.mechanical.bulk_modulus',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.mechanical.bulk_modulus.type',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.mechanical.bulk_modulus.value',
                                    ),
                                ],
                            ),
                            MenuItemNestedObject(
                                path='results.properties.mechanical.shear_modulus',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.mechanical.shear_modulus.type',
                                        show_input=False,
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.mechanical.shear_modulus.value',
                                    ),
                                ],
                            ),
                            MenuItemNestedObject(
                                path='results.properties.mechanical.energy_volume_curve',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.mechanical.energy_volume_curve.type',
                                        options=5,
                                    ),
                                ],
                            ),
                        ]
                    ),
                    'usecases': Menu(),
                    'solarcell': Menu(
                        items=[
                            MenuItemHistogram(
                                x='results.properties.optoelectronic.solar_cell.efficiency',
                            ),
                            MenuItemHistogram(
                                x='results.properties.optoelectronic.solar_cell.fill_factor',
                            ),
                            MenuItemHistogram(
                                x='results.properties.optoelectronic.solar_cell.open_circuit_voltage',
                            ),
                            MenuItemHistogram(
                                x='results.properties.optoelectronic.solar_cell.short_circuit_current_density',
                            ),
                            MenuItemHistogram(
                                x='results.properties.optoelectronic.solar_cell.illumination_intensity',
                            ),
                            MenuItemHistogram(
                                x='results.properties.optoelectronic.solar_cell.device_area',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.device_architecture',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.device_stack',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.absorber',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.absorber_fabrication',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.electron_transport_layer',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.hole_transport_layer',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.substrate',
                            ),
                            MenuItemTerms(
                                search_quantity='results.properties.optoelectronic.solar_cell.back_contact',
                            ),
                        ]
                    ),
                    'heterogeneouscatalyst': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.properties.catalytic.reaction.name',
                            ),
                            MenuItemNestedObject(
                                path='results.properties.catalytic.reaction.reactants',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.catalytic.reaction.reactants.name',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.catalytic.reaction.reactants.conversion',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.catalytic.reaction.reactants.gas_concentration_in',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.catalytic.reaction.reactants.gas_concentration_out',
                                    ),
                                ],
                            ),
                            MenuItemNestedObject(
                                path='results.properties.catalytic.reaction.products',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.catalytic.reaction.products.name',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.catalytic.reaction.products.selectivity',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.catalytic.reaction.products.gas_concentration_out',
                                    ),
                                ],
                            ),
                            MenuItemHistogram(
                                x='results.properties.catalytic.reaction.reaction_conditions.temperature',
                            ),
                            MenuItemNestedObject(
                                path='results.properties.catalytic.catalyst',
                                items=[
                                    MenuItemTerms(
                                        search_quantity='results.properties.catalytic.catalyst.catalyst_type',
                                    ),
                                    MenuItemTerms(
                                        search_quantity='results.properties.catalytic.catalyst.preparation_method',
                                    ),
                                    MenuItemTerms(
                                        search_quantity='results.properties.catalytic.catalyst.catalyst_name',
                                    ),
                                    MenuItemTerms(
                                        search_quantity='results.properties.catalytic.catalyst.characterization_methods',
                                    ),
                                    MenuItemHistogram(
                                        x='results.properties.catalytic.catalyst.surface_area',
                                    ),
                                ],
                            ),
                        ]
                    ),
                    'author': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='authors.name',
                                options=0,
                            ),
                            MenuItemHistogram(
                                x='upload_create_time',
                            ),
                            MenuItemTerms(
                                search_quantity='external_db',
                                options=5,
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='datasets.dataset_name',
                            ),
                            MenuItemTerms(
                                search_quantity='datasets.doi',
                                options=0,
                            ),
                        ]
                    ),
                    'metadata': Menu(
                        items=[
                            MenuItemVisibility(),
                            MenuItemTerms(
                                search_quantity='entry_id',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='upload_id',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='upload_name',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.material.material_id',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='datasets.dataset_id',
                                options=0,
                            ),
                            MenuItemDefinitions(),
                        ]
                    ),
                    'optimade': Menu(items=[MenuItemOptimade()]),
                    'eln': Menu(
                        items=[
                            MenuItemTerms(
                                search_quantity='results.eln.sections',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.eln.tags',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.eln.methods',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.eln.instruments',
                                show_input=False,
                            ),
                            MenuItemTerms(
                                search_quantity='results.eln.names',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.eln.descriptions',
                                options=0,
                            ),
                            MenuItemTerms(
                                search_quantity='results.eln.lab_ids',
                                options=0,
                            ),
                        ]
                    ),
                    'custom_quantities': Menu(items=[MenuItemCustomQuantities()]),
                    'combine': MenuItemTerms(
                        search_quantity='combine',
                        options={
                            True: MenuItemOption(
                                label='Combine results from several entries',
                                description='If selected, your filters may be matched from several entries that contain the same material. When unchecked, the material has to have a single entry that matches all your filters.',
                            )
                        },
                        show_header=False,
                        show_input=False,
                        show_statistics=False,
                    ),
                }.get(key)
                if not menu:
                    continue
                size = value.get('size')
                new_size = {
                    FilterMenuSizeEnum.S: MenuSizeEnum.MD,
                    FilterMenuSizeEnum.M: MenuSizeEnum.LG,
                    FilterMenuSizeEnum.L: MenuSizeEnum.XL,
                    FilterMenuSizeEnum.XL: MenuSizeEnum.XXL,
                    None: MenuSizeEnum.MD,
                }.get(size)
                if isinstance(menu, Menu):
                    menu.title = value.get('label')
                    if new_size:
                        menu.size = new_size
                    menu.indentation = value.get('level')
                items.append(menu)
            del values['filter_menus']
            values['menu'] = Menu(title='Filters', size=MenuSizeEnum.SM, items=items)

        return values


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
