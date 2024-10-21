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
        columns=[
            Column(quantity=f'data.name#{schema}', selected=True),
            Column(quantity=f'data.population#{schema}', selected=True),
            Column(quantity=f'data.area#{schema}', selected=True),
            Column(quantity=f'data.population_density#{schema}', selected=True),
            Column(quantity=f'data.literacy#{schema}', selected=True),
            Column(quantity=f'data.net_migration#{schema}', selected=True),
            Column(quantity=f'data.infant_mortality#{schema}', selected=True),
            Column(quantity=f'data.birthrate#{schema}', selected=True),
            Column(quantity=f'data.deathrate#{schema}', selected=True),
        ],
        filters=Filters(include=[f'*{schema}']),
        filter_menus=FilterMenus(options={}),
        dashboard=Dashboard(
            widgets=[
                WidgetScatterPlot(
                    layout={'lg': Layout(h=6, w=8, x=0, y=0)},
                    x=Axis(quantity=f'data.literacy#{schema}'),
                    y=Axis(quantity=f'data.industry#{schema}'),
                    size='1000',
                    autorange=True,
                ),
                WidgetScatterPlot(
                    layout={'lg': Layout(h=6, w=8, x=8, y=0)},
                    x=Axis(quantity=f'data.literacy#{schema}'),
                    y=Axis(quantity=f'data.agriculture#{schema}'),
                    size='1000',
                    autorange=True,
                ),
                WidgetScatterPlot(
                    layout={'lg': Layout(h=6, w=8, x=16, y=0)},
                    x=Axis(quantity=f'data.literacy#{schema}'),
                    y=Axis(quantity=f'data.service#{schema}'),
                    size='1000',
                    autorange=True,
                ),
                WidgetHistogram(
                    layout={'lg': Layout(h=3, w=8, x=0, y=6)},
                    quantity=f'data.phones#{schema}',
                    scale='1/2',
                    nbins=30,
                    show_input=False,
                    autorange=True,
                ),
                WidgetHistogram(
                    layout={'lg': Layout(h=3, w=8, x=8, y=6)},
                    quantity=f'data.birthrate#{schema}',
                    scale='linear',
                    nbins=30,
                    show_input=False,
                    autorange=True,
                ),
                WidgetHistogram(
                    layout={'lg': Layout(h=3, w=8, x=16, y=6)},
                    quantity=f'data.net_migration#{schema}',
                    scale='linear',
                    nbins=30,
                    show_input=False,
                    autorange=True,
                ),
            ]
        ),
    ),
)
