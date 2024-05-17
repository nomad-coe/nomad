from nomad.config.models.plugins import AppEntryPoint
from nomad.config.models.ui import (
    App,
    Axis,
    Column,
    Columns,
    Dashboard,
    FilterMenus,
    Filters,
    Layout,
    WidgetHistogram,
    WidgetScatterPlot,
)

schema = 'nomad_countries.schema_packages.country.Country'

country = AppEntryPoint(
    name='Countries of the World app',
    app=App(
        label='Countries of the World',
        path='countries',
        description='Search app for the tutorial data',
        category='Other',
        filters_locked={'section_defs.definition_qualified_name': schema},
        columns=Columns(
            selected=[
                f'data.name#{schema}',
                f'data.population#{schema}',
                f'data.area#{schema}',
                f'data.population_density#{schema}',
                f'data.literacy#{schema}',
                f'data.net_migration#{schema}',
                f'data.infant_mortality#{schema}',
                f'data.birthrate#{schema}',
                f'data.deathrate#{schema}',
            ],
            options={
                f'data.name#{schema}': Column(),
                f'data.population#{schema}': Column(),
                f'data.area#{schema}': Column(),
                f'data.population_density#{schema}': Column(),
                f'data.literacy#{schema}': Column(),
                f'data.net_migration#{schema}': Column(),
                f'data.infant_mortality#{schema}': Column(),
                f'data.birthrate#{schema}': Column(),
                f'data.deathrate#{schema}': Column(),
            },
        ),
        filters=Filters(include=[f'*{schema}']),
        filter_menus=FilterMenus(options={}),
        dashboard=Dashboard(
            widgets=[
                WidgetScatterPlot(
                    type='scatterplot',
                    layout={'lg': Layout(h=6, w=8, x=0, y=0)},
                    x=Axis(quantity=f'data.literacy#{schema}'),
                    y=Axis(quantity=f'data.industry#{schema}'),
                    size='1000',
                    autorange=True,
                ),
                WidgetScatterPlot(
                    type='scatterplot',
                    layout={'lg': Layout(h=6, w=8, x=8, y=0)},
                    x=Axis(quantity=f'data.literacy#{schema}'),
                    y=Axis(quantity=f'data.agriculture#{schema}'),
                    size='1000',
                    autorange=True,
                ),
                WidgetScatterPlot(
                    type='scatterplot',
                    layout={'lg': Layout(h=6, w=8, x=16, y=0)},
                    x=Axis(quantity=f'data.literacy#{schema}'),
                    y=Axis(quantity=f'data.service#{schema}'),
                    size='1000',
                    autorange=True,
                ),
                WidgetHistogram(
                    type='histogram',
                    layout={'lg': Layout(h=3, w=8, x=0, y=6)},
                    quantity=f'data.phones#{schema}',
                    scale='1/2',
                    nbins=30,
                    showinput=False,
                    autorange=True,
                ),
                WidgetHistogram(
                    type='histogram',
                    layout={'lg': Layout(h=3, w=8, x=8, y=6)},
                    quantity=f'data.birthrate#{schema}',
                    scale='linear',
                    nbins=30,
                    showinput=False,
                    autorange=True,
                ),
                WidgetHistogram(
                    type='histogram',
                    layout={'lg': Layout(h=3, w=8, x=16, y=6)},
                    quantity=f'data.net_migration#{schema}',
                    scale='linear',
                    nbins=30,
                    showinput=False,
                    autorange=True,
                ),
            ]
        ),
    ),
)
