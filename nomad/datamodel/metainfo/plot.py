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
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import Quantity, SubSection, Package, MSection, JSON
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np
from copy import deepcopy
from datetime import datetime


class PlotlyError(Exception):
    """Plotly errors."""

    pass


class PlotSectionError(Exception):
    """PlotSection errors."""

    pass


def get_figure_layout(annotation):
    label = annotation.get('label', None)
    if label is not None:
        annotation.pop('label')
    index = annotation.get('index', None)
    if index is not None:
        annotation.pop('index')
    return label, index


def express_do_plot(plotly_express_annotation, archive, logger):
    method_name = plotly_express_annotation.pop('method')
    layout = plotly_express_annotation.get('layout', None)
    if layout:
        plotly_express_annotation.pop('layout')
    traces = plotly_express_annotation.get('traces', [])
    if traces:
        plotly_express_annotation.pop('traces')
    method = getattr(px, method_name)
    kwargs = {}
    for key, value in plotly_express_annotation.items():
        if isinstance(value, list):
            items = []
            for item in value:
                if item.startswith('#'):
                    resolved_value = None
                    try:
                        resolved_value = archive.m_resolve(f'/data/{item[1:]}')
                        if (
                            resolved_value
                            and resolved_value[0]
                            and isinstance(resolved_value[0], datetime)
                        ):
                            resolved_value = [
                                value.isoformat() for value in resolved_value
                            ]
                    except Exception as e:
                        logger.warning(e)
                    items.append(resolved_value)
                else:
                    items.append(item)
            kwargs[key] = items
        elif isinstance(value, str):
            if value.startswith('#'):
                resolved_value = None
                try:
                    resolved_value = archive.m_resolve(f'/data/{value[1:]}')
                    if (
                        resolved_value
                        and resolved_value[0]
                        and isinstance(resolved_value[0], datetime)
                    ):
                        resolved_value = [value.isoformat() for value in resolved_value]
                except Exception as e:
                    logger.warning(e)
                kwargs[key] = resolved_value
            else:
                kwargs[key] = value
    try:
        figure = method(**kwargs)
        return figure, layout, traces
    except Exception as e:
        raise PlotlyError(e)


def convert_to_list(data):
    x = data.get('x', None)
    y = data.get('y', None)
    z = data.get('z', None)
    marker = data.get('marker', None)
    if x is not None and isinstance(x, np.ndarray):
        data['x'] = x.tolist()
    if y is not None and isinstance(y, np.ndarray):
        data['y'] = y.tolist()
    if z is not None and isinstance(z, np.ndarray):
        data['z'] = z.tolist()
    if marker is not None:
        color = marker.get('color', None)
        if color is not None and isinstance(color, np.ndarray):
            marker['color'] = color.tolist()


m_package = Package()


class Figure(MSection):
    label = Quantity(type=str, description='Label shown in the plot selection.')
    index = Quantity(type=int, description='Index of figure in the plot selection.')


class PlotlyFigureQuantity(Quantity):
    def __set__(self, obj, value):
        if obj is None:
            raise KeyError(
                'Cannot overwrite quantity definition. Only values can be set.'
            )

        # Make generated json serializable by converting numpy.ndarray to python list
        if value is not None:
            if 'data' in value:
                all_data = value['data']
                if isinstance(all_data, list):
                    for data in all_data:
                        convert_to_list(data)
                else:
                    convert_to_list(all_data)

        obj.m_set(self, value)


class PlotlyFigure(Figure):
    figure = PlotlyFigureQuantity(
        type=JSON, description='Contains the JSON serialization for a plotly figure.'
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.figure = kwargs.get('figure', None)


class PlotSection(ArchiveSection):
    """
    The PlotSection base section serves as an additional functionality to your sections.
    This base section is designed to simplify the process of creating various types of
    plots, making it easy to use Plotly Express, Plotly Subplot, and the general Plotly graph objects.

    Features:
    - Plotly Express: Create simple and quick plots with a high-level, expressive API.
    - Plotly Subplot: Organize multiple plots into subplots for more complex visualizations.
    - General Plotly Graph Objects: Fine-tune your plots by working directly with Plotly's graph objects.

    Usage:
    - Inherit from this base section to leverage its plot functionality.
    - Customize your plots using the annotations plotly-express, plotly-subplots, or/and plotly-graph-object.

    Example Usage:
    ```python
    class CustomSection(PlotSection, ElnBaseSection, EntryData):
        m_def = Section()
        time = Quantity(type=float, shape=['*'], unit='s', a_eln=dict(component='NumberEditQuantity'))
        substrate_temperature = Quantity(type=float, shape=['*'], unit='K', a_eln=dict(component='NumberEditQuantity'))
        chamber_pressure = Quantity(type=float, shape=['*'], unit='Pa', a_eln=dict(component='NumberEditQuantity'))

        def normalize(self, archive, logger):
            super(CustomSection, self).normalize(archive, logger)

            first_line = px.scatter(x=self.time, y=self.substrate_temperature)
            second_line = px.scatter(x=self.time, y=self.chamber_pressure)
            figure1 = make_subplots(rows=1, cols=2, shared_yaxes=True)
            figure1.add_trace(first_line.data[0], row=1, col=1)
            figure1.add_trace(second_line.data[0], row=1, col=2)
            figure1.update_layout(height=400, width=716, title_text="Creating Subplots in Plotly")
            self.figures.append(PlotlyFigure(label='figure 1', figure=figure1.to_plotly_json()))

            figure2 = px.scatter(x=self.substrate_temperature, y=self.chamber_pressure, color=self.chamber_pressure, title="Chamber as a function of Temperature")
            self.figures.append(PlotlyFigure(label='figure 2', index=1, figure=figure2.to_plotly_json()))

            heatmap_data = [[None, None, None, 12, 13, 14, 15, 16],
                 [None, 1, None, 11, None, None, None, 17],
                 [None, 2, 6, 7, None, None, None, 18],
                 [None, 3, None, 8, None, None, None, 19],
                 [5, 4, 10, 9, None, None, None, 20],
                 [None, None, None, 27, None, None, None, 21],
                 [None, None, None, 26, 25, 24, 23, 22]]

            heatmap = go.Heatmap(z=heatmap_data, showscale=False, connectgaps=True, zsmooth='best')
            figure3 = go.Figure(data=heatmap)
            self.figures.append(PlotlyFigure(label='figure 3', index=0, figure=figure3.to_plotly_json()))
    ```
    """

    figures = SubSection(
        sub_section=PlotlyFigure.m_def, repeats=True, label_quantity='label'
    )

    def normalize(self, archive, logger):
        super(PlotSection, self).normalize(archive, logger)

        all_figures = []
        plotly_express_annotations = deepcopy(
            self.m_def.m_get_annotations('plotly_express', None)
        )
        plotly_graph_object_annotations = deepcopy(
            self.m_def.m_get_annotations('plotly_graph_object', None)
        )
        plotly_subplots_annotations = deepcopy(
            self.m_def.m_get_annotations('plotly_subplots', None)
        )

        if isinstance(plotly_express_annotations, dict):
            plotly_express_annotations = [plotly_express_annotations]
        if isinstance(plotly_graph_object_annotations, dict):
            plotly_graph_object_annotations = [plotly_graph_object_annotations]
        if isinstance(plotly_subplots_annotations, dict):
            plotly_subplots_annotations = [plotly_subplots_annotations]

        if plotly_graph_object_annotations:
            for plotly_graph_object_annotation in plotly_graph_object_annotations:
                figure = {}
                label, figure_index = get_figure_layout(plotly_graph_object_annotation)
                if 'data' in plotly_graph_object_annotation:
                    figure['data'] = plotly_graph_object_annotation['data']
                if 'layout' in plotly_graph_object_annotation:
                    figure['layout'] = plotly_graph_object_annotation['layout']
                if 'config' in plotly_graph_object_annotation:
                    figure['config'] = plotly_graph_object_annotation['config']
                all_figures.append(
                    {'label': label, 'index': figure_index, 'graph_object': figure}
                )

        if plotly_express_annotations:
            for plotly_express_annotation in plotly_express_annotations:
                try:
                    label, figure_index = get_figure_layout(plotly_express_annotation)
                    fig, layout, traces = express_do_plot(
                        plotly_express_annotation, archive, logger
                    )
                    try:
                        all_traces = go.Figure(fig.data[0])
                        total_layout = layout
                    except Exception as e:
                        raise PlotlyError(e)
                    for trace in traces:
                        trace_fig, layout, _ = express_do_plot(trace, archive, logger)
                        try:
                            all_traces.add_trace(trace_fig.data[0])
                        except Exception as e:
                            raise PlotlyError(e)
                    if total_layout:
                        try:
                            all_traces.update_layout(**total_layout)
                        except Exception as e:
                            raise PlotlyError(e)
                    plotly_graph_object = all_traces.to_plotly_json()
                    all_figures.append(
                        {
                            'label': label,
                            'index': figure_index,
                            'graph_object': plotly_graph_object,
                        }
                    )
                except Exception as error:
                    if isinstance(error, PlotlyError):
                        # Handle internal plotly errors as warning to not stop all figures. Easier to figure out problem.
                        logger.warning(error)
                    else:
                        raise PlotSectionError(error)

        if plotly_subplots_annotations:
            for plotly_subplots_annotation in plotly_subplots_annotations:
                try:
                    label, figure_index = get_figure_layout(plotly_subplots_annotation)
                    parameters = plotly_subplots_annotation.get('parameters', None)
                    if parameters:
                        rows = parameters.get('rows', None)
                        cols = parameters.get('cols', None)
                        if rows and cols:
                            subplot_titles = [''] * rows * cols
                            plotly_express_list = plotly_subplots_annotation.get(
                                'plotly_express', []
                            )
                            for index, plotly_express in enumerate(plotly_express_list):
                                title = plotly_express.get('title', None)
                                if title is not None:
                                    title = (
                                        plotly_express.get('layout', {})
                                        .get('title', {})
                                        .get('text', None)
                                    )
                                if title is not None:
                                    subplot_titles[index] = title

                            default_parameters = {
                                **{'subplot_titles': subplot_titles},
                                **parameters,
                            }
                            try:
                                figure = make_subplots(**default_parameters)
                            except Exception as e:
                                raise PlotlyError(e)

                            for row in range(1, rows + 1):
                                for col in range(1, cols + 1):
                                    if len(plotly_express_list) > 0:
                                        plotly_express_annotation = (
                                            plotly_express_list.pop(0)
                                        )
                                        fig, sub_layout, traces = express_do_plot(
                                            plotly_express_annotation, archive, logger
                                        )
                                        try:
                                            if sub_layout:
                                                fig.update_layout(**sub_layout)
                                            figure.add_trace(
                                                fig.data[0], row=row, col=col
                                            )
                                            for trace in traces:
                                                trace_fig, layout, _ = express_do_plot(
                                                    trace, archive, logger
                                                )
                                                try:
                                                    figure.add_trace(
                                                        trace_fig.data[0],
                                                        row=row,
                                                        col=col,
                                                    )
                                                except Exception as e:
                                                    raise PlotlyError(e)
                                            if sub_layout:
                                                xaxis = sub_layout.get('xaxis', None)
                                                if xaxis:
                                                    figure.update_xaxes(
                                                        xaxis, row=row, col=col
                                                    )
                                                yaxis = sub_layout.get('yaxis', None)
                                                if yaxis:
                                                    figure.update_yaxes(
                                                        yaxis, row=row, col=col
                                                    )
                                        except Exception as e:
                                            raise PlotlyError(e)
                            layout = plotly_subplots_annotation.get('layout', None)
                            if layout:
                                try:
                                    figure.update_layout(**layout)
                                except Exception as e:
                                    raise PlotlyError(e)
                            plotly_graph_object = figure.to_plotly_json()
                            all_figures.append(
                                {
                                    'label': label,
                                    'index': figure_index,
                                    'graph_object': plotly_graph_object,
                                }
                            )
                except Exception as error:
                    if isinstance(error, PlotlyError):
                        # Handle internal plotly errors as warning to not stop all figures. Easier to figure out problem.
                        logger.warning(error)
                    else:
                        raise PlotSectionError(error)

        # Only remove existing figure, if we add more figures.
        # Ideally, we would just add figures. But with our current ELN implementation
        # this might add the same set of figures over and over again.
        # But, if we always remove all existing figures, figures produced and uploaded
        # by users, will be lost. This is a bit of a tradeoff.
        if len(all_figures) > 0:
            self.figures = []

        for figure in all_figures:
            self.figures.append(
                PlotlyFigure(
                    label=figure['label'],
                    index=figure['index'],
                    figure=figure['graph_object'],
                )
            )


m_package.__init_metainfo__()
