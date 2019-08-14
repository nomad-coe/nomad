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

from .client import client


def codes(client):
    result = client.repo.search(per_page=1, owner='admin').response().result

    return sorted([
        code for code, values in result.quantities['code_name'].items()
        if code != 'not processed' and values['code_runs'] > 0], key=lambda x: x.lower())


def error_fig(client):
    labels = codes(client)

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


def codes_fig(client):
    labels = codes(client)

    # get the data
    result = client.repo.search(
        per_page=1, owner='admin', metrics=['total_energies']).response().result
    all_entries = {
        code: values['code_runs']
        for code, values in result.quantities['code_name'].items()
        if code != 'not processed' and code in labels}
    total_energies = {
        code: values['total_energies']
        for code, values in result.quantities['code_name'].items()
        if code != 'not processed' and code in labels}

    fig, ax1 = plt.subplots(figsize=(15, 6), dpi=72)
    x = np.arange(len(labels))  # the label locations
    width = 0.7 / 2  # the width of the bars
    plt.sca(ax1)
    plt.xticks(rotation=90)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_title('Number of entries (code runs, sets of input/output files) and total energy calculations per code')

    color = 'tab:red'
    ax1.set_yscale('power', exponent=0.25)
    ax1.set_ylabel('number of entries', color=color)
    ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.1f}M'))
    ax1.bar(x - width / 2, [all_entries[code] / 1e6 for code in labels], width, label='entries', color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_yscale('power', exponent=0.25)
    ax2.set_ylabel('number of total energy calculations', color=color)
    ax2.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.1f}M'))
    ax2.bar(x + width / 2, [total_energies[code] / 1e6 for code in labels], width, label='total energy calculations', color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()
    plt.show()


@client.command(help='Generate various matplotlib charts')
@click.option('--errors', is_flag=True, help='Two charts with relative and absolute parser/normalizer errors per code.')
@click.option('--per-code', is_flag=True, help='Entries and total energy calculations per code.')
def statistics(errors, per_code):
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

    if errors:
        error_fig(client)
    if per_code:
        codes_fig(client)
