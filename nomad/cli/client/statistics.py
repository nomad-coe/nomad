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

"""
A command that generates various statistics.
"""

from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import click
import json
from datetime import datetime

from .client import client


def codes(client, minimum=1, **kwargs):
    data = client.repo.search(per_page=1, **kwargs).response().result

    x_values = sorted([
        code for code, values in data.statistics['code_name'].items()
        if code != 'not processed' and values.get('calculations', 1000) >= minimum], key=lambda x: x.lower())

    return data.statistics, x_values, 'code_name', 'code'


def dates(client, minimum=1, **kwargs):
    data = client.repo.search(per_page=1, date_histogram=True, **kwargs).response().result

    x_values = list([
        x for x in data.statistics['date_histogram'].keys()])

    return data.statistics, x_values, 'date_histogram', 'month'


def error_fig(client):
    _, labels, _, _ = codes(client)

    def code_values(metric='code_runs', **kwargs):
        result = client.repo.search(
            per_page=1,
            owner='admin',
            metrics=[] if metric == 'code_runs' else metric,
            **kwargs).response().result

        return {
            code: values[metric]
            for code, values in result.quantities['code_name'].items()
            if code != 'not processed' and (not labels or code in labels) > 0}

    # get the data
    all_entries = code_values()
    parser_failure_label = 'parser failure'
    error_types = [
        {'name': parser_failure_label, 'search': dict(system='not processed')},
        {'name': 'failed system classification', 'search': dict(system='unavailable')},
        {'name': 'no basis set available', 'search': dict(basis_set='unavailable')},
        {'name': 'no XC functional available', 'search': dict(xc_functional='unavailable')}
    ]
    errors = {
        error_type['name']: {
            code: failures
            for code, failures in code_values(**error_type['search']).items()}
        for error_type in error_types}
    errors_rates = {
        error_type['name']: {
            code: 0 if all_entries[code] == 0 else failures / all_entries[code]
            for code, failures in code_values(**error_type['search']).items()}
        for error_type in error_types}

    fig, axs = plt.subplots(figsize=(15, 12), dpi=72, nrows=2)

    def draw_error_chart(errors, ax, colors, entries=None, mul=1, scale=0.5):
        n_bars = len(errors) - 1
        leg_colors = list(colors)

        x = np.arange(len(labels))  # the label locations
        width = 0.7 / n_bars  # the width of the bars
        plt.sca(ax)
        plt.xticks(rotation=90)

        if entries is not None:
            ax.bar(x, [entries[code] for code in labels], width * n_bars, label='all entries', color=colors.pop(0))

        i = -1
        not_processed = [errors[parser_failure_label][code] * mul for code in labels]
        ax.bar(x, not_processed, width * n_bars, label=parser_failure_label, color=colors.pop(0))
        for key, values in errors.items():
            if key != parser_failure_label:
                ax.bar(x + i * width, [values[code] * mul for code in labels], width, label=key, bottom=not_processed, color=colors.pop(0))
                i += 1

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_yscale('power', exponent=scale)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        leg = ax.get_legend()
        for i in range(0, len(leg_colors)):
            leg.legendHandles[i].set_color(leg_colors[i])

        fig.tight_layout()

    ax = axs[0]
    ax.set_title('Absolute number of entries with parser errors or missing repository metadata compared to all entries per code')
    ax.set_ylabel('number of entries', )
    colors = ['grey', 'red', 'yellow', 'orange', 'brown']
    draw_error_chart(errors, ax, entries=all_entries, mul=1, scale=0.25, colors=colors)
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

    ax = axs[1]
    ax.set_title('Relative rates of entries with parser errors or missing repository metadata per code')
    ax.set_ylabel('rate in %', )
    colors = ['red', 'yellow', 'orange', 'brown']
    draw_error_chart(errors_rates, ax, mul=100, colors=colors)

    plt.show()

    return fig, plt


class Metric:
    def __init__(self, metric, label=None, power=None, multiplier=1, format=None, cumulate=False):
        if label is None:
            label = metric

        self.metric = metric
        self.agg = None
        self.label = label
        self.multiplier = multiplier
        self.power = power
        self.format = format
        self.cumulate = cumulate

    def draw_axis(self, axis, data, x_values, x_positions, width, color, only=False):
        label_color = 'black' if only else color
        value_map = {
            x: values[self.metric]
            for x, values in data[self.agg].items()
            if x in x_values}

        if self.power is not None:
            axis.set_yscale('power', exponent=self.power)
        else:
            axis.set_yscale('log')
        axis.set_ylabel(self.label, color=label_color)

        if self.format is not None:
            axis.yaxis.set_major_formatter(ticker.StrMethodFormatter(self.format))

        y_values = [value_map[x] * self.multiplier for x in x_values]
        if self.cumulate:
            y_values = np.array(y_values).cumsum()
        axis.bar(x_positions, y_values, width, label=self.label, color=color, align='edge')
        axis.tick_params(axis='y', labelcolor=label_color)

        # TODO remove
        # for x, v in zip(x_positions, y_values):
        #     axis.text(x + .1, v, ' {:,}'.format(int(v)), color=color, fontweight='bold', rotation=90)
        # import matplotlib.lines as mlines
        # line = mlines.Line2D([min(x_positions), max(x_positions)], [72, 72], color=color)
        # axis.add_line(line)


def bar_plot(
        client, retrieve, metric1, metric2=None, title=None, format_xlabel=None,
        xlim={}, ylim=dict(bottom=1), **kwargs):
    if format_xlabel is None:
        format_xlabel = lambda x: x

    metrics = [] if metric1.metric == 'code_runs' else [metric1.metric]
    if metric2 is not None:
        metrics += [] if metric2.metric == 'code_runs' else [metric2.metric]

    data, x_values, agg, agg_label = retrieve(client, metrics=metrics, statistics=True, **kwargs)
    metric1.agg = agg
    if metric2 is not None:
        metric2.agg = agg

    fig, ax1 = plt.subplots(figsize=(5, 4), dpi=72)
    x = np.arange(len(x_values))
    width = 0.8 / 2
    if metric2 is None:
        width = 0.8
    plt.sca(ax1)
    plt.xticks(rotation=90)
    ax1.set_xticks(x)
    ax1.set_xticklabels([format_xlabel(value) if value != 'Quantum Espresso' else 'Q. Espresso' for value in x_values])
    ax1.margins(x=0.01)
    ax1.set_xlim(**xlim)
    # i = 0
    # for label in ax1.xaxis.get_ticklabels():
    #     label.set_visible(i % 4 == 0)
    #     i += 1

    if title is None:
        title = 'Number of %s' % metric1.label
        if metric2 is not None:
            title += ' and %s' % metric2.label
        title += ' per %s' % agg_label
        ax1.set_title(title)
    elif title != '':
        ax1.set_title(title)

    metric1.draw_axis(ax1, data, x_values, x - (width / 2), width, 'tab:blue', only=metric2 is None)
    ax1.set_ylim(**ylim)
    ax1.set_yticks([40, 30, 20, 10, 5, 1, 0.5, 0.1])
    ax1.grid(which='major', axis='y', linestyle='--')

    if metric2:
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        metric2.draw_axis(
            ax2, data, x_values, x + width / 2, width, 'tab:red')
        ax2.set_ylim(bottom=1)

    fig.tight_layout()

    return fig, plt


@client.command(help='Generate various matplotlib charts')
@click.option('--errors', is_flag=True, help='Two charts with relative and absolute parser/normalizer errors per code.')
@click.option('--x-axis', type=str, help='Aggregation used for x-axis, values are "code" and "time".')
@click.option('--y-axis', multiple=True, type=str, help='Metrics used for y-axis, values are "entries", "energies", "calculations", "users".')
@click.option('--cumulate', is_flag=True, help='Cumulate over x-axis.')
@click.option('--title', type=str, help='Override chart title with given value.')
@click.option('--total', is_flag=True, help='Provide total sums of key metrics.')
@click.option('--save', type=str, help='Save to given file instead of showing the plot.')
@click.option('--power', type=float, help='User power scale instead of log with the given inverse power.')
@click.option('--open-access', is_flag=True, help='Only consider Open-Access data.')
@click.option('--minimum', type=int, default=1, help='Only consider codes with at least the given ammount of entries.')
def statistics(errors, title, x_axis, y_axis, cumulate, total, save, power, open_access, minimum):
    from .client import create_client
    client = create_client()

    class PowerScale(mscale.ScaleBase):
        name = 'power'

        def __init__(self, axis, exponent, **kwargs):
            mscale.ScaleBase.__init__(self, axis, **kwargs)
            self.exponent = exponent

        def set_default_locators_and_formatters(self, axis):
            axis.set_major_locator(ticker.AutoLocator())
            axis.set_major_formatter(ticker.ScalarFormatter())
            axis.set_minor_locator(ticker.NullLocator())
            axis.set_minor_formatter(ticker.NullFormatter())

        def limit_range_for_scale(self, vmin, vmax, minpos):
            return max(0., vmin), vmax

        class Transform(mtransforms.Transform):
            input_dims = 1
            output_dims = 1
            is_separable = True

            def __init__(self, exponent):
                super().__init__()
                self.exponent = exponent

            def transform_non_affine(self, a):
                return np.array(a)**self.exponent

            def inverted(self):
                return PowerScale.Transform(1 / self.exponent)

        def get_transform(self):
            return self.Transform(self.exponent)

    mscale.register_scale(PowerScale)

    kwargs = {}
    if cumulate:
        kwargs.update(
            power=1,
            multiplier=1e-6,
            format='{x:,.1f}M')
    elif power is not None:
        kwargs.update(
            power=1 / power,
            multiplier=1e-6,
            format='{x:,.1f}M')

    metrics = {
        'entries': Metric(
            'code_runs',
            label='entries (code runs)',
            cumulate=cumulate,
            **kwargs),
        'users': Metric(
            'users',
            cumulate=cumulate,
            label='users that provided data'),
        'energies': Metric(
            'total_energies',
            label='total energy calculations',
            cumulate=cumulate,
            **kwargs),
        'calculations': Metric(
            'calculations',
            label='calculations (e.g. total energy)',
            cumulate=cumulate,
            **kwargs)
    }

    if errors:
        fig, plt = error_fig(client)

    owner = 'all' if open_access else 'admin'

    if x_axis is not None:
        assert 1 <= len(y_axis) <= 2, 'Need 1 or 2 y axis'

        kwargs = {}
        if x_axis == 'code':
            x_axis = codes
            kwargs.update(ylim=dict(bottom=0))
        elif x_axis == 'time':
            x_axis = dates
            kwargs.update(
                ylim=dict(bottom=0),
                format_xlabel=lambda x: datetime.fromtimestamp(int(x) / 1000).strftime('%b %y'))
        else:
            assert False, 'x axis can only be "code" or "time"'

        y_axis = [metrics[y] for y in y_axis]

        fig, plt = bar_plot(
            client, x_axis, *y_axis, title=title, owner=owner, minimum=minimum, **kwargs)

    if errors or x_axis is not None:
        if save is not None:
            fig.savefig(save, bbox_inches='tight')
        else:
            plt.show()

    if total:
        data = client.repo.search(
            per_page=1, owner=owner, statistics=True,
            metrics=['total_energies', 'calculations', 'uploaders', 'authors', 'datasets']).response().result
        print(json.dumps(data.statistics['total'], indent=4))
