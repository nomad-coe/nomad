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

import numpy as np
import ase

from nomad.metainfo import Package, Quantity, SubSection, Datetime, MEnum, SectionProxy
from nomad.datamodel.results import Results, Material
from nomad.datamodel.data import EntryData


m_package = Package(name='eln_examples')


class ElnBaseSection(EntryData):
    name = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'))
    description = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity', props=dict(multiline=True, minRows=10)))


class Experiment(ElnBaseSection):
    pass


class Sample(ElnBaseSection):
    sample_id = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'))
    chemical_formula = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'))

    sample_type = Quantity(
        type=MEnum('slap', 'bulk', 'library'))

    weight = Quantity(
        type=np.dtype(np.float64), unit='kg',
        a_eln=dict(component='NumberEditQuantity'))
    grain_size = Quantity(
        type=np.dtype(np.float64), unit='m',
        a_eln=dict(component='NumberEditQuantity'))
    cycles = Quantity(
        type=np.dtype(np.int32),
        a_eln=dict(component='NumberEditQuantity'))

    coating_method = Quantity(
        type=MEnum(
            'Vapor deposition', 'Chemical vapor deposition', 'Metalorganic vapour phase epitaxy',
            'Electrostatic spray assisted vapour deposition (ESAVD)', 'Sherardizing',
            'Some forms of Epitaxy', 'Molecular beam epitaxy', 'Physical vapor deposition',
            'Cathodic arc deposition', 'Electron beam physical vapor deposition (EBPVD)',
            'Ion plating', 'Ion beam assisted deposition (IBAD)', 'Magnetron sputtering',
            'Pulsed laser deposition', 'Sputter deposition', 'Vacuum deposition',
            'Pulsed electron deposition (PED)', 'Chemical and electrochemical techniques',
            'Conversion coating', 'Autophoretic', 'Anodising',
            'Chromate conversion coating', 'Plasma electrolytic oxidation', 'Phosphate',
            'Ion beam mixing', 'Pickled and oiled, a type of plate steel coating',
            'Plating', 'Electroless plating', 'Electroplating', 'Spraying', 'Spray painting',
            'High velocity oxygen fuel (HVOF)', 'Plasma spraying',
            'Thermal spraying', 'Kinetic metallization (KM)', 'Plasma transferred wire arc thermal spraying',
            'The common forms of Powder coating', 'Roll-to-roll coating processes',
            'Common roll-to-roll coating processes include:', 'Air knife coating', 'Anilox coater',
            'Flexo coater', 'Gap Coating', 'Knife-over-roll coating', 'Gravure coating',
            'Immersion dip coating', 'Kiss coating', 'Metering rod (Meyer bar) coating',
            'Roller coating', 'Forward roller coating', 'Reverse roll coating',
            'Silk Screen coater', 'Rotary screen', 'Lithography', 'Flexography',
            'Physical coating processes', 'Langmuir-Blodgett', 'Spin coating', 'Dip coating'),
        a_eln=dict(
            label='Coating Method',
            component='AutocompleteEditQuantity'))

    sample_status = Quantity(
        type=MEnum('prepared', 'characterized', 'depleted', 'unknown'),
        a_eln=dict(
            label='Sample Status',
            component='EnumEditQuantity'))

    spin_polarized = Quantity(
        type=bool,
        a_eln=dict(
            label='Spin Polarized',
            component='BoolEditQuantity'))

    processes = SubSection(section_def=SectionProxy('Process'), repeats=True)
    experiment = SubSection(section_def=SectionProxy('Experiment'), repeats=False)

    def normalize_results(self, results: Results, logger):
        if self.chemical_formula:
            if not results.material:
                results.material = Material()
            try:
                atoms = ase.Atoms(self.chemical_formula)
                results.material.chemical_formula_hill = atoms.get_chemical_formula(mode='hill')
                results.material.chemical_formula_reduced = atoms.get_chemical_formula(mode='reduce')
                results.material.chemical_formula_descriptive = results.material.chemical_formula_hill
                results.material.elements = list(set(atoms.get_chemical_symbols()))
            except Exception as e:
                logger.warn('could not normalize formula', exc_info=e)


class Process(ElnBaseSection):
    start_time = Quantity(
        type=Datetime,
        a_eln=dict(component='StringEditQuantity'))


m_package.__init_metainfo__()
