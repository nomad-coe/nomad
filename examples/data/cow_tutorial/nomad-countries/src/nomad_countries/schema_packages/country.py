import numpy as np
import plotly.express as px
import pandas as pd
from nomad.metainfo import MSection, Quantity, SubSection, Package
from nomad.datamodel import Schema
from nomad.datamodel.metainfo.plot import PlotSection, PlotlyFigure

m_package = Package(name='Countries of the World')


class Timeseries(MSection):
    year = Quantity(type=np.int, shape=['*'])
    value = Quantity(type=np.float64, shape=['*'])


class Country(PlotSection, Schema):
    name = Quantity(type=str)
    population = Quantity(type=np.int32)
    area = Quantity(type=np.float64, unit='km^2')
    population_density = Quantity(type=np.float64, unit='1/km^2')
    coastline = Quantity(type=np.float64, description='cost/area ratio')
    net_migration = Quantity(type=np.float64)
    infant_mortality = Quantity(type=np.float64, description='per 1000 births')
    literacy = Quantity(
        type=np.float64, description='Literacy in % of adult population'
    )
    phones = Quantity(type=np.float64, description='Phones per 1,000 people')
    birthrate = Quantity(type=np.float64, description='per 1,000 people per year')
    deathrate = Quantity(type=np.float64, description='per 1,000 people per year')
    agriculture = Quantity(type=np.float64)
    industry = Quantity(type=np.float64)
    service = Quantity(type=np.float64)

    gdp = SubSection(
        section=Timeseries, description='GDP per capita (constant 2005 US$)'
    )
    birth_rate = SubSection(section=Timeseries, description='per 1,000 people per year')

    def normalize(self, archive, logger):
        super(Country, self).normalize(archive, logger)
        archive.metadata.entry_name = self.name
        self.population_density = self.population / self.area

        self.figures.append(
            PlotlyFigure(
                figure=px.line(
                    pd.DataFrame(dict(year=self.gdp.year, GDP=self.gdp.value)),
                    x='year',
                    y=['GDP'],
                ).to_plotly_json()
            )
        )


m_package.__init_metainfo__()
