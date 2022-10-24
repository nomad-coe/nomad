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
from nomad import utils
from nomad.units import ureg
from nomad.datamodel.data import EntryData, ArchiveSection, user_reference, author_reference
from nomad.metainfo.metainfo import SectionProxy
from nomad.datamodel.results import ELN, Results, Material, BandGap
from nomad.metainfo import Package, Quantity, Datetime, Reference, Section
from nomad.datamodel.metainfo.eln.perovskite_solar_cell_database import addSolarCell

m_package = Package(name='material_library')


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
        a_eln=dict(component='StringEditQuantity'))

    lab_id = Quantity(
        type=str,
        description='A id string that is unique at least for the lab that produced this data.',
        a_eln=dict(component='StringEditQuantity'))

    description = Quantity(
        type=str,
        description=(
            'A humand description. This provides room for human readable information '
            'that could not be captured in the ELN.'),
        a_eln=dict(component='RichTextEditQuantity'))

    def normalize(self, archive, logger):
        super(ElnBaseSection, self).normalize(archive, logger)

        if isinstance(self, EntryData):
            if archive.data == self and self.name:
                archive.metadata.entry_name = self.name
            EntryData.normalize(self, archive, logger)

        if not archive.results:
            archive.results = Results(eln=ELN())
        if not archive.results.eln:
            archive.results.eln = ELN()

        if self.lab_id:
            if archive.results.eln.lab_ids is None:
                archive.results.eln.lab_ids = []
            archive.results.eln.lab_ids.append(self.lab_id)

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


class ElnActivityBaseSection(ElnBaseSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used for
    laboratory activities, e.g. processes, characterizations, measurements, etc.
    '''
    datetime = Quantity(
        type=Datetime,
        description='The date and time when this activity was done.',
        a_eln=dict(component='DateTimeEditQuantity'))

    method = Quantity(
        type=str,
        description='A short consistent handle for the applied method.')

    user = Quantity(
        type=author_reference,
        description='The corresponding user for the activity.',
        a_eln=dict(component='AuthorEditQuantity'))

    def normalize(self, archive, logger):
        super(ElnActivityBaseSection, self).normalize(archive, logger)

        if archive.results.eln.methods is None:
            archive.results.eln.methods = []
        if self.method:
            archive.results.eln.methods.append(self.method)
        else:
            archive.results.eln.methods.append(self.m_def.name)


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
        from ase import Atoms
        from pymatgen.core import Composition
        if self.chemical_formula:
            if not archive.results:
                archive.results = Results()
            if not archive.results.material:
                archive.results.material = Material()
            material = archive.results.material

            try:
                pycom = Composition(self.chemical_formula).get_integer_formula_and_factor()[0]
                atoms = Atoms(pycom)
                material.elements = list(set(atoms.get_chemical_symbols()))
                material.chemical_formula_hill = atoms.get_chemical_formula(mode='hill')
                material.chemical_formula_reduced = atoms.get_chemical_formula(mode='reduce')
                material.chemical_formula_descriptive = self.chemical_formula
            except Exception as e:
                logger.warn('could not analyse chemical formula', exc_info=e)


class Chemical(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for chemicals.'''
    pass


class Sample(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for samples.'''
    pass


class Instrument(ElnBaseSection):
    ''' A ELN base section that can be used for instruments.'''
    def normalize(self, archive, logger):
        super(Instrument, self).normalize(archive, logger)

        if self.name:
            if archive.results.eln.instruments is None:
                archive.results.eln.instruments = []
            archive.results.eln.instruments.append(self.name)


class Process(ElnActivityBaseSection):
    ''' A ELN base section that can be used for processes.'''
    pass


class Measurement(ElnActivityBaseSection):
    ''' A ELN base section that can be used for measurements.'''
    pass


class SampleID(ArchiveSection):
    ''' A ELN base section that can be used for sample IDs.
    If the `sample_owner`, `sample_short_name`, `ìnstitute`, and `creation_datetime`
    quantities are provided, the sample_id will be automatically created as a combination
    of these four quantities.'''

    # m_def = Section(
    #     a_eln=dict(hide=['name', 'lab_id']))

    institute = Quantity(
        type=str,
        description='Alias/short name of the home institute of the owner, i.e. *HZB*.',
        a_eln=dict(component='StringEditQuantity'))

    sample_owner = Quantity(
        type=str,
        shape=[],
        description='Name or alias of the process operator, e.g. jmp',
        a_eln=dict(component='StringEditQuantity'))

    creation_datetime = Quantity(
        type=Datetime,
        description='Creation date of the sample.',
        a_eln=dict(component='DateTimeEditQuantity'))

    sample_short_name = Quantity(
        type=str,
        description='''A short name of the sample (the identifier scribed on the smaple,
         or in the sample container), e.g. 4001-8, YAG-2-34.
         This is to be managed and decided internally by the labs,
         although we recomend to avoid the following characters on it: "_", "/", "\" and "."''',
        a_eln=dict(component='StringEditQuantity'))

    sample_id = Quantity(
        type=str,
        description='''Full sample id. Ideally a human readable sample id convention,
        which is simple, understandable and still having chances of becoming unique.
        If the `sample_owner`, `sample_short_name`, `ìnstitute`, and `creation_datetime`
        are provided, this will be formed automatically by joining these components by an underscore (_).
        Spaces in any of the individual components will be replaced with hyphens (-).
        An example would be hzb_oah_20200602_4001-08''',
        a_eln=dict(component='StringEditQuantity'))

    children = Quantity(
        type=Reference(SectionProxy('SampleID')),
        shape=["*"],
        descriptions='A reference to a sample which are children of this one.',
        a_eln=dict(component='ReferenceEditQuantity'))

    parents = Quantity(
        type=Reference(SectionProxy('SampleID')),
        descriptions='A reference to sample which are parents of this one.',
        a_eln=dict(component='ReferenceEditQuantity'))

    def normalize(self, archive, logger):
        super(SampleID, self).normalize(archive, logger)

        if self.institute and self.sample_short_name and self.sample_owner and self.creation_datetime:
            creation_date = self.creation_datetime.strftime('%Y%m%d')
            sample_owner = self.sample_owner.replace(' ', '-')
            sample_id_list = [self.institute, sample_owner, creation_date, self.sample_short_name]
            self.sample_id = '_'.join(sample_id_list)

        if isinstance(self, EntryData):
            if archive.data == self and self.sample_id:
                archive.metadata.entry_name = self.sample_id.replace('_', ' ')
            EntryData.normalize(self, archive, logger)

        if not archive.results:
            archive.results = Results(eln=ELN())
        if not archive.results.eln:
            archive.results.eln = ELN()

        if self.sample_id:
            if not archive.results.eln.lab_ids:
                archive.results.eln.lab_ids = []
            archive.results.eln.lab_ids.append(self.sample_id.replace('_', ' '))

        if not archive.results.eln.sections:
            archive.results.eln.sections = []
        archive.results.eln.sections.append(self.m_def.name)


class PublicationReference(ArchiveSection):
    ''' A ELN base section that can be used for references.'''

    DOI_number = Quantity(
        type=str,
        shape=[],
        description="""
            The DOI number referring to the published paper or dataset where the data can be found.
            Examples:
            10.1021/jp5126624
            10.1016/j.electacta.2017.06.032
        """,
        a_eln=dict(
            component='EnumEditQuantity', props=dict(suggestions=[])))

    publication_authors = Quantity(
        type=str,
        shape=['*'],
        description="""
            The authors of the publication.
            If several authors, end with et al. If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """)

    publication_date = Quantity(
        type=Datetime,
        shape=[],
        description="""
            Publication date.
            If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """)

    journal = Quantity(
        type=str,
        shape=[],
        description="""
            Name of the journal where the data is published.
            If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """)

    publication_title = Quantity(
        type=str,
        shape=[],
        description="""
            Title of the publication.
            If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """)

    def normalize(self, archive, logger):
        super(PublicationReference, self).normalize(archive, logger)
        from nomad.datamodel.datamodel import EntryMetadata
        import dateutil.parser
        import requests

        # Parse journal name, lead author and publication date from crossref
        if self.DOI_number:
            try:
                url = f'https://api.crossref.org/works/{self.DOI_number}?mailto=contact@nomad-lab.eu'
                timeout = 2
                r = requests.get(url, timeout=timeout)
                if r.status_code == 200:
                    temp_dict = r.json()
                    self.publication_authors = [f"{v['given']} {v['family']}" for v in temp_dict['message']['author']]
                    self.journal = temp_dict['message']['container-title'][0]
                    self.publication_title = temp_dict['message']['title'][0]
                    self.publication_date = dateutil.parser.parse(temp_dict['message']['created']['date-time'])
                    if not archive.metadata:
                        archive.metadata = EntryMetadata()
                    if not archive.metadata.references:
                        archive.metadata.references = []
                    if self.DOI_number not in archive.metadata.references:
                        archive.metadata.references.append(self.DOI_number)

                else:
                    logger.warning(f'Could not parse DOI number {self.DOI_number}')
            except Exception as e:
                logger.warning(f'Could not parse crossref for {self.DOI_number}')
                logger.warning(e)


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
        addSolarCell(archive)
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
        addSolarCell(archive)
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
        addSolarCell(archive)
        if self.bandgap is not None:
            band_gap = BandGap()
            band_gap.value = self.bandgap
            if band_gap.value is None:
                band_gap.value = 0
            archive.results.properties.optoelectronic.band_gap = [band_gap]
            props = archive.results.properties.available_properties
            if not props:
                props = []
            props.append('optoelectronic.band_gap')
            archive.results.properties.available_properties = props


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

        addSolarCell(archive)
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
                self.urbach_energy = eqe_dict['urbach_e']
                self.photon_energy_array = np.array(eqe_dict['interpolated_photon_energy'])
                self.raw_photon_energy_array = np.array(eqe_dict['photon_energy_raw'])
                self.eqe_array = np.array(eqe_dict['interpolated_eqe'])
                self.raw_eqe_array = np.array(eqe_dict['eqe_raw'])

        if self.photon_energy_array is not None:
            self.wavelength_array = self.photon_energy_array.to('nm', 'sp')  # pylint: disable=E1101
            self.raw_wavelength_array = self.raw_photon_energy_array.to('nm', 'sp')  # pylint: disable=E1101

        addSolarCell(archive)
        if self.bandgap_eqe is not None:
            band_gap = BandGap()
            band_gap.value = self.bandgap_eqe
            if band_gap.value is None:
                band_gap.value = 0
            archive.results.properties.optoelectronic.band_gap = [band_gap]
            props = archive.results.properties.available_properties
            if not props:
                props = []
            props.append('optoelectronic.band_gap')
            archive.results.properties.available_properties = props


m_package.__init_metainfo__()
