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
import datetime
import re
from typing import Any, Dict, List
from nomad import utils
from nomad.units import ureg
from nomad.datamodel.data import EntryData, ArchiveSection, author_reference, BasicElnCategory
from nomad.metainfo.metainfo import MSection, MProxy, MEnum, Category, MCategory
from nomad.datamodel.results import ELN, Results, Material
from nomad.datamodel.results import ElementalComposition as ResultsElementalComposition
from nomad.metainfo import Package, Quantity, Datetime, Reference, Section, SubSection
from ase.data import chemical_symbols, atomic_numbers, atomic_masses
from nomad.datamodel.metainfo.eln.perovskite_solar_cell_database import (
    add_solar_cell, add_band_gap
)

from nomad.datamodel.metainfo.simulation.system import System as SystemTheory, Atoms
from nomad.datamodel.metainfo.simulation.run import Run
from nomad.datamodel.metainfo.eln.nexus_data_converter import (
    NexusDataConverter,
    ElnYamlConverter
)

from nomad.datamodel.metainfo.basesections import (
    Entity,
    Activity,
    ElementalComposition,
    System,
    Instrument,
    Process,
    Measurement,
    SynthesisMethod,
    Component,
    Ensemble,
    CASExperimentalProperty,
    CASPropertyCitation,
    CASSubstance as Substance,
    SampleID,
    PublicationReference,
)


m_package = Package(name='eln')


class User(MSection):
    user = Quantity(
        type=author_reference,
        description='The corresponding user for the activity.',
        a_eln=dict(component='AuthorEditQuantity'))


class ElnBaseSection(ArchiveSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used properties.

    If you inherit from this section, but do not need some quantities, list those
    quantities in the `eln.hide` annotation of your inheriting section definition.

    Besides predefining some quantities, these base sections will add some metadata
    to NOMAD's search. A particular example are `tags`, if you define a string
    or enum quantity in your sections named `tags`, its values will be searchable.
    '''

    name = Quantity(
        type=str,
        description='A short human readable and descriptive name.',
        a_eln=dict(component='StringEditQuantity', label='Short name'))

    datetime = Quantity(
        type=Datetime,
        description='The date and time associated with this section.',
        a_eln=dict(component='DateTimeEditQuantity'))

    lab_id = Quantity(
        type=str,
        description='''An ID string that is unique at least for the lab that produced this
            data.''',
        a_eln=dict(component='StringEditQuantity', label="ID"))

    description = Quantity(
        type=str,
        description='Any information that cannot be captured in the other fields.',
        a_eln=dict(component='RichTextEditQuantity'))

    def normalize(self, archive, logger):
        super(ElnBaseSection, self).normalize(archive, logger)

        if isinstance(self, EntryData):
            if archive.data == self and self.name:
                archive.metadata.entry_name = self.name
            elif self.name is None:
                self.name = archive.metadata.entry_name.split('.')[0].replace('_', ' ')
            EntryData.normalize(self, archive, logger)

        if not archive.results:
            archive.results = Results(eln=ELN())
        if not archive.results.eln:
            archive.results.eln = ELN()

        if self.datetime is None:
            self.datetime = datetime.datetime.now()

        if self.lab_id:
            if archive.results.eln.lab_ids is None:
                archive.results.eln.lab_ids = []
            if self.lab_id not in archive.results.eln.lab_ids:
                archive.results.eln.lab_ids.append(self.lab_id)
        else:
            if archive.results.eln.lab_ids and len(archive.results.eln.lab_ids) > 0:
                self.lab_id = archive.results.eln.lab_ids[0]

        if getattr(self, 'name'):
            if archive.results.eln.names is None:
                archive.results.eln.names = []
            archive.results.eln.names.append(self.name)

        if getattr(self, 'description'):
            if archive.results.eln.descriptions is None:
                archive.results.eln.descriptions = []
            archive.results.eln.descriptions.append(self.description)

        if getattr(self, 'tags', None):
            if archive.results.eln.tags is None:
                archive.results.eln.tags = []
            tags = self.tags
            if isinstance(tags, list):
                archive.results.eln.tags.extend(tags)
            else:
                archive.results.eln.tags.append(tags)

        if not archive.results.eln.sections:
            archive.results.eln.sections = []
        archive.results.eln.sections.append(self.m_def.name)


class BasicEln(ElnBaseSection, EntryData):
    ''' The most basic ELN to instantiate. '''
    m_def = Section(categories=[BasicElnCategory], a_eln=dict(lane_width='600px'), label='Basic ELN')

    tags = Quantity(
        type=str,
        shape=['*'],
        description='Add a tag that can be used for search.',
        a_eln=dict(component='StringEditQuantity'))


# Legacy sections:
class ElnWithFormulaBaseSection(ElnBaseSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used for
    items with a chemical formula, e.g. chemicals or samples.
    '''
    chemical_formula = Quantity(
        type=str,
        description=(
            'The chemical formula. This will be used directly and '
            'indirectly in the search. The formula will be used itself as well as '
            'the extracted chemical elements.'),
        a_eln=dict(component='StringEditQuantity'))

    def normalize(self, archive, logger):
        super(ElnWithFormulaBaseSection, self).normalize(archive, logger)

        if logger is None:
            logger = utils.get_logger(__name__)
        from nomad.atomutils import Formula
        if self.chemical_formula:
            if not archive.results:
                archive.results = Results()
            if not archive.results.material:
                archive.results.material = Material()
            try:
                formula = Formula(self.chemical_formula)
            except Exception as e:
                logger.warn('could not analyse chemical formula', exc_info=e)
            else:
                try:
                    formula.populate(archive.results.material)
                except ValueError as e:
                    logger.info('composition information already defined, skipping populating it based on formula', exc_info=e)


class Chemical(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for chemicals.'''
    pass


class Sample(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for samples.'''
    pass


class ElnWithStructureFile(ArchiveSection):
    '''
    A base section for for parsing crystal structure files, e.g. `.cif`, and
    populating the Material section in Results.
    '''

    structure_file = Quantity(
        type=str,
        description='The structure file.',
        a_eln=dict(component='FileEditQuantity'))

    def normalize(self, archive, logger):
        super(ElnWithStructureFile, self).normalize(archive, logger)

        if self.structure_file:
            from ase.io import read
            from nomad.normalizing.results import ResultsNormalizer
            from nomad.normalizing.system import SystemNormalizer
            from nomad.normalizing.optimade import OptimadeNormalizer
            from nomad.datamodel.metainfo.simulation.run import Program

            with archive.m_context.raw_file(self.structure_file) as f:
                try:
                    structure = read(f.name)
                except Exception as e:
                    raise ValueError('could not read structure file') from e

                run = Run()
                archive.run = [run]
                system = SystemTheory()
                system.atoms = Atoms()
                try:
                    system.atoms.lattice_vectors = structure.get_cell() * ureg.angstrom
                except Exception as e:
                    logger.warn('Could not parse lattice vectors from structure file.', exc_info=e)
                system.atoms.labels = structure.get_chemical_symbols()
                system.atoms.positions = structure.get_positions() * ureg.angstrom
                try:
                    system.atoms.periodic = structure.get_pbc()
                except Exception as e:
                    logger.warn('Could not parse periodicity from structure file.', exc_info=e)
                    system.atoms.periodic = [True, True, True]
                system.atoms.species = structure.get_atomic_numbers()
                archive.run[0].system = [system]
                program = Program()
                archive.run[0].program = program
                archive.run[0].program.name = 'Structure File Importer'
                system_normalizer = SystemNormalizer(archive)
                system_normalizer.normalize()
                optimade_normalizer = OptimadeNormalizer(archive)
                optimade_normalizer.normalize()
                results_normalizer = ResultsNormalizer(archive)
                results_normalizer.normalize()

        # TODO: rewrite it in a way in which the run section is not needed and System is
        # directly added to the archive.data
        # set run to None if exist
        if archive.run:
            archive.run = None


# Solar cell sections:
class SolarCellDefinition(ArchiveSection):

    stack_sequence = Quantity(
        type=str,
        shape=['*'],
        description="""
            The stack sequence describing the cell. Use the following formatting guidelines
            - Start with the substrate to the left and list the materials in each layer of the device
            - If two materials, e.g. A and B, are mixed in one layer, list the materials in alphabetic order and separate them with semicolons, as in (A; B)
            - The absorber layer in other databases is commonly stated with a generaic name as “Perovskite”, regardless of composition, mixtures, dimensionality etc.
                There are other fields to describe in depth the absorber layer.
        """,
        a_eln=dict(
            component='EnumEditQuantity', props=dict(suggestions=sorted([]))))

    solar_cell_area = Quantity(
        type=np.dtype(np.float64),
        unit='cm**2',
        shape=[],
        description="""
            The total cell area in cm^2.
            The total area is defined as the area that would provide photovoltaic performance.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    architecture = Quantity(
        type=str,
        shape=[],
        description="""
            The cell architecture with respect to the direction of current flow and
            the order in which layers are deposited.
            The two most common are nip (also referred to as normal) and pin (also referred to as inverted)
            but there are also a few others, e.g. Back contacted.
            - *nip* architecture means that the electrons are collected at the substrate side.
            The typical example is in perovskite solar cells when a TiO2 electron selective contact is deposited
            between the perovskite and the substrate (e.g. SLG | FTO | TiO2-c | Perovskite | …)
            - *pin* architecture means that it instead is the holes that are collected at the substrate side. The typical example is when a PEDOT:PSS hole selective contact is deposited between the perovskite and the substrate (e.g. SLG | FTO | PEDOT:PSS |Perovskite | …)
        """,
        a_eln=dict(
            component='EnumEditQuantity',
            props=dict(
                suggestions=['Unknown', 'Pn-Heterojunction', 'Front contacted', 'Back contacted', 'pin', 'nip', 'Schottky'])))

    def normalize(self, archive, logger):
        super(SolarCellDefinition, self).normalize(archive, logger)
        add_solar_cell(archive)
        if self.stack_sequence:
            if '/' in self.stack_sequence:
                archive.results.properties.optoelectronic.solar_cell.device_stack = self.stack_sequence.split('/')
            elif '|' in self.stack_sequence:
                archive.results.properties.optoelectronic.solar_cell.device_stack = self.stack_sequence.split(' | ')
            else:
                archive.results.properties.optoelectronic.solar_cell.device_stack = self.stack_sequence

        if self.architecture:
            archive.results.properties.optoelectronic.solar_cell.device_architecture = self.architecture
        if self.solar_cell_area:
            archive.results.properties.optoelectronic.solar_cell.device_area = self.solar_cell_area
        if not archive.results.material:
            archive.results.material = Material()
        material = archive.results.material
        if material.functional_type is None:
            material.functional_type = ['semiconductor', 'solar cell']


class SolarCellLayer(ArchiveSection):

    solar_cell_layer_type = Quantity(
        type=str,
        shape=[],
        description='type of the layer',
        a_eln=dict(component='EnumEditQuantity',
                   props=dict(suggestions=[
                    'Substrate',
                    'Absorber',
                    'Hole Transport Layer',
                    'Electron Transport Layer',
                    'Contact',
                    'Buffer',
                    'p-type contact',
                    'n-type contact',
                    'other'])))

    layer_name = Quantity(
        type=str,
        shape=['0..*'],
        description="""
            The name of the layer.
        """,
        a_eln=dict(
            component='EnumEditQuantity',
            props=dict(suggestions=[])))

    layer_thickness = Quantity(
        type=np.dtype(np.float64),
        unit=('nm'),
        shape=[],
        description="""
            The thickness of the layer in nm.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    def normalize(self, archive, logger):
        super(SolarCellLayer, self).normalize(archive, logger)
        add_solar_cell(archive)
        if self.layer_name:
            if self.solar_cell_layer_type == 'Absorber':
                archive.results.properties.optoelectronic.solar_cell.absorber = [self.layer_name]
            elif self.solar_cell_layer_type == 'Substrate':
                archive.results.properties.optoelectronic.solar_cell.substrate = [self.layer_name]
            elif self.solar_cell_layer_type == 'Hole Transport Layer':
                archive.results.properties.optoelectronic.solar_cell.hole_transport_layer = [self.layer_name]
            elif self.solar_cell_layer_type == 'Electron Transport Layer':
                archive.results.properties.optoelectronic.solar_cell.electron_transport_layer = [self.layer_name]
            elif self.solar_cell_layer_type == 'Contact':
                archive.results.properties.optoelectronic.solar_cell.back_contact = [self.layer_name]


class SolarCellBaseSectionWithOptoelectronicProperties(ArchiveSection):

    bandgap = Quantity(
        type=np.dtype(np.float64),
        unit=('eV'),
        shape=[],
        description="""
            The bandgap of the solar cell.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    def normalize(self, archive, logger):
        super(SolarCellBaseSectionWithOptoelectronicProperties, self).normalize(archive, logger)
        add_solar_cell(archive)
        add_band_gap(archive, self.bandgap)


class SolarCellJV(ArchiveSection):

    m_def = Section(label_quantity='cell_name')

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    certified_values = Quantity(
        type=bool,
        shape=[],
        description="""
            TRUE if the IV data is measured by an independent and certification institute.
            If your solar simulator is calibrated by a calibrated reference diode,
            that does not count as a certified result.
        """,
        a_eln=dict(
            component='BoolEditQuantity'))

    certification_institute = Quantity(
        type=str,
        shape=[],
        description="""
            The name of the certification institute that has measured the certified device.
            Example:
            Newport
            NIM, National Institute of Metrology of China
            KIER, Korea Institute of Energy Research
        """,
        a_eln=dict(
            component='EnumEditQuantity', props=dict(suggestions=sorted([
                'National Institute ofMetrology, China',
                'Quality supervision＆Testing Center of Chemical＆Physical Power Sources of Information Industry',
                'CREST, Photovoltaic Meaasurement and calibration Laboratory at Universit of Loughborough',
                'Photovoltaic and Wind Power Systems Quality Test Center, Chinese Academy of Sciences',
                'NREL', 'Institute of Metrology (NIM) of China',
                'PVEVL, National Central University, Taiwan',
                'NIM, National Institute of Metrology of China',
                'Fraunhofer ISE',
                'SIMIT, Shanghai Institute of Microsystem and Information Technology',
                'Newport',
                'CSIRO, PV Performance Lab at Monash University',
                'AIST, National Institute of Advanced Industrial Science and Technology',
                'CPVT, National Center of Supervision and Inspection on Solar Photovoltaic Products Quality of China',
                'KIER, Korea Institute of Energy Research',
                'Newport Corporation',
                'Solar Power Lab at Arizona State University']))))

    light_intensity = Quantity(
        type=np.dtype(np.float64),
        unit=('mW/cm**2'),
        shape=[],
        default=100.0,
        description="""
            The light intensity during the IV measurement
            - If there are uncertainties, only state the best estimate, e.g. write 100 and not 90-100.
            - Standard AM 1.5 illumination correspond to 100 mW/cm2
            - If you need to convert from illumination given in lux; at 550 nm, 1 mW/cm2 corresponds to 6830 lux. Be aware that the conversion change with the spectrum used. As a rule of thumb for general fluorescent/LED light sources, around 0.31mW corresponded to 1000 lux. If your light intensity is measured in lux, it probably means that your light spectra deviates quite a lot from AM 1.5, wherefore it is very important that you also specify the light spectra in the next column.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    open_circuit_voltage = Quantity(
        type=np.dtype(np.float64),
        unit='V',
        shape=[],
        description="""
            Open circuit voltage.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    short_circuit_current_density = Quantity(
        type=np.dtype(np.float64),
        unit='mA / cm**2',
        shape=[],
        description="""
            Short circuit current density.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    fill_factor = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description="""
            Fill factor.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    efficiency = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description="""
            Power conversion efficiency.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    potential_at_maximum_power_point = Quantity(
        type=np.dtype(np.float64),
        unit='V',
        shape=[],
        description="""
            The potential at the maximum power point, Vmp.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    current_density_at_maximun_power_point = Quantity(
        type=np.dtype(np.float64),
        unit='mA / cm**2',
        shape=[],
        description="""
            The current density at the maximum power point, *Jmp*.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    series_resistance = Quantity(
        type=np.dtype(np.float64),
        unit='ohm*cm**2',
        shape=[],
        description="""
            The series resistance as extracted from the *J-V* curve.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    shunt_resistance = Quantity(
        type=np.dtype(np.float64),
        unit='ohm*cm**2',
        shape=[],
        description="""
            The shunt resistance as extracted from the *J-V* curve.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    def derive_n_values(self):
        if self.current_density is not None:
            return len(self.current_density)
        if self.voltage is not None:
            return len(self.voltage)
        else:
            return 0

    n_values = Quantity(type=int, derived=derive_n_values)

    def update_results(self, archive):
        if self.open_circuit_voltage is not None:
            archive.results.properties.optoelectronic.solar_cell.open_circuit_voltage = self.open_circuit_voltage
        if self.short_circuit_current_density is not None:
            archive.results.properties.optoelectronic.solar_cell.short_circuit_current_density = self.short_circuit_current_density
        if self.fill_factor is not None:
            archive.results.properties.optoelectronic.solar_cell.fill_factor = self.fill_factor
        if self.efficiency is not None:
            archive.results.properties.optoelectronic.solar_cell.efficiency = self.efficiency
        if self.light_intensity is not None:
            archive.results.properties.optoelectronic.solar_cell.illumination_intensity = self.light_intensity

    def normalize(self, archive, logger):
        super(SolarCellJV, self).normalize(archive, logger)

        add_solar_cell(archive)
        self.update_results(archive)


class SolarCellJVCurve(SolarCellJV):

    def cell_params(self):
        """
        Calculates basic solar cell parametes form a current density (mA/cm**2)
        voltage (V) curve.

        Returns:
            Voc (V) open circuit voltage
            Jsc (mA/cm**2) short circuit current density
            FF fill factor in absolute values (0-1)
            efficiency power conversion efficiency in percentage (0-100)
        """
        from scipy import interpolate
        j_v_interpolated = interpolate.interp1d(self.current_density, self.voltage)
        Voc = j_v_interpolated(0)
        v_j_interpolated = interpolate.interp1d(self.voltage, self.current_density)
        Isc = v_j_interpolated(0)
        if Isc >= 0:
            idx = np.argmax(self.voltage * self.current_density)
        else:
            idx = np.argmin(self.voltage * self.current_density)
        Vmp = self.voltage[idx]
        Imp = self.current_density[idx]
        Isc = abs(Isc)
        FF = abs(Vmp.magnitude * Imp.magnitude / (Voc * Isc))
        efficiency = Voc * FF * Isc
        return Voc, Isc, FF, efficiency

    cell_name = Quantity(
        type=str,
        shape=[],
        description='Cell identification name.',
        a_eln=dict(component='StringEditQuantity'))

    current_density = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='mA/cm^2',
        description='Current density array of the *JV* curve.',
        a_plot={
            'x': 'voltage', 'y': 'current_density'
        })

    voltage = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='V',
        description='Voltage array of the of the *JV* curve.',
        a_plot={
            'x': 'voltage', 'y': 'current_density'
        })

    def normalize(self, archive, logger):
        super(SolarCellJVCurve, self).normalize(archive, logger)
        if self.current_density is not None:
            if self.voltage is not None:
                self.open_circuit_voltage, self.short_circuit_current_density, self.fill_factor, self.efficiency = self.cell_params()
                self.update_results(archive)


class SolarCellEQE(ArchiveSection):

    m_def = Section(
        a_eln=dict(lane_width='600px'))

    eqe_data_file = Quantity(
        type=str,
        description="""
                    Drop here your eqe file and click save for processing.
                    """,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    header_lines = Quantity(
        type=np.dtype(np.int64),
        default=0,
        description="""
        Number of header lines in the file. Edit in case your file has a header.
        """,
        a_eln=dict(component='NumberEditQuantity'))

    light_bias = Quantity(
        type=np.dtype(np.float64),
        unit=('mW/cm**2'),
        shape=[],
        description="""
        The light intensity of any bias light during the EQE measurement.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    bandgap_eqe = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='eV',
        description="""
        Bandgap derived from the EQE spectrum.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    integrated_jsc = Quantity(
        type=np.dtype(np.float64),
        unit='mA / cm**2',
        shape=[],
        description="""
        The integrated short circuit current density $J_{SC}$ from the product of the EQE spectrum
        with the *AM 1.5G* sun spectrum.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    integrated_j0rad = Quantity(
        type=np.dtype(np.float64),
        unit='mA / cm**2',
        shape=[],
        description="""
        The integrated $J_{0, Rad}$ derived from the EQE data.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    voc_rad = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='V',
        description="""
        Radiative $V_{OC}$ derived from the EQE data in V.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    urbach_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='eV',
        description="""
        Urbach energy fitted from the eqe in eV.
        """,
        a_eln=dict(
            component='NumberEditQuantity'))

    urbach_energy_fit_std_dev = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='eV',
        description="""
        Standard deviation of the fitted Urbach energy parameter from the eqe in eV.
        """)

    def derive_n_values(self):
        if self.eqe_array is not None:
            return len(self.eqe_array)
        if self.photon_energy_array is not None:
            return len(self.photon_energy_array)
        else:
            return 0

    n_values = Quantity(type=int, derived=derive_n_values)

    def derive_n_raw_values(self):
        if self.raw_eqe_array is not None:
            return len(self.raw_eqe_array)
        if self.raw_photon_energy_array is not None:
            return len(self.raw_photon_energy_array)
        else:
            return 0

    n_raw_values = Quantity(type=int, derived=derive_n_raw_values)

    raw_eqe_array = Quantity(
        type=np.dtype(np.float64), shape=['n_raw_values'],
        description='EQE array of the spectrum',
        a_plot={
            'x': 'photon_energy_array', 'y': 'raw_eqe_array'
        })

    raw_photon_energy_array = Quantity(
        type=np.dtype(np.float64), shape=['n_raw_values'], unit='eV',
        description='Raw Photon energy array of the eqe spectrum',
        a_plot={
            'x': 'raw_photon_energy_array', 'y': 'raw_eqe_array'
        })

    raw_wavelength_array = Quantity(
        type=np.dtype(np.float64), shape=['n_raw_values'], unit='nanometer',
        description='Raw wavelength array of the eqe spectrum',
        a_plot={
            'x': 'raw_wavelength_array', 'y': 'raw_eqe_array'
        })

    eqe_array = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        description='EQE array of the spectrum',
        a_plot={
            'x': 'photon_energy_array', 'y': 'eqe_array'
        })

    wavelength_array = Quantity(
        type=np.dtype(np.float64), shape=['n_values'], unit='nanometer',
        description='Interpolated/extrapolated wavelength array with *E<sub>u</sub>* of the eqe spectrum ',
        a_plot={
            'x': 'wavelength_array', 'y': 'eqe_array'
        })

    photon_energy_array = Quantity(
        type=np.dtype(np.float64), shape=['n_values'], unit='eV',
        description='Interpolated/extrapolated photon energy array with a *E<sub>u</sub>*  of the eqe spectrum',
        a_plot={
            'x': 'photon_energy_array', 'y': 'eqe_array'
        })

    def normalize(self, archive, logger):
        super(SolarCellEQE, self).normalize(archive, logger)
        from .perovskite_solar_cell_database.eqe_parser import EQEAnalyzer
        if (self.eqe_data_file):
            with archive.m_context.raw_file(self.eqe_data_file) as f:
                eqe_dict = EQEAnalyzer(f.name, header_lines=self.header_lines).eqe_dict()
                self.measured = True
                self.bandgap_eqe = eqe_dict['bandgap']
                self.integrated_jsc = eqe_dict['jsc'] * ureg('A/m**2')
                self.integrated_j0rad = eqe_dict['j0rad'] * ureg('A/m**2') if 'j0rad' in eqe_dict else logger.warning('The j0rad could not be calculated.')
                self.voc_rad = eqe_dict['voc_rad'] if 'voc_rad' in eqe_dict else logger.warning('The voc_rad could not be calculated.')
                self.urbach_energy = eqe_dict['urbach_e'] if 'urbach_e' in eqe_dict else logger.warning('The urbach_energy could not be calculated.')
                if 'error_urbach_std' in eqe_dict:
                    self.urbach_energy_fit_std_dev = eqe_dict['error_urbach_std']
                self.photon_energy_array = np.array(eqe_dict['interpolated_photon_energy'])
                self.raw_photon_energy_array = np.array(eqe_dict['photon_energy_raw'])
                self.eqe_array = np.array(eqe_dict['interpolated_eqe'])
                self.raw_eqe_array = np.array(eqe_dict['eqe_raw'])

        if self.photon_energy_array is not None:
            self.wavelength_array = self.photon_energy_array.to('nm', 'sp')  # pylint: disable=E1101
            self.raw_wavelength_array = self.raw_photon_energy_array.to('nm', 'sp')  # pylint: disable=E1101

        add_solar_cell(archive)
        add_band_gap(archive, self.bandgap_eqe)


m_package.__init_metainfo__()
