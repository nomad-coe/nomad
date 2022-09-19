import numpy as np
from nomad.metainfo import MSection, Quantity, SubSection


class Element(MSection):
    label = Quantity(type=str)
    density = Quantity(type=np.float64, unit='g/cm**3')
    isotopes = Quantity(type=np.int32, shape=['*'])


class Sample(MSection):
    composition = Quantity(type=str)
    elements = SubSection(section=Element, repeats=True)
