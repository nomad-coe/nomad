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

from nomad.units import ureg
from nomad.metainfo import (
    MSection, Package, Quantity, SubSection, MEnum, Reference, Datetime, Section)
from nomad.datamodel.data import EntryData, ArchiveSection, ElnExampleCategory


m_package = Package(name='material_library')


class Chemical(EntryData):
    '''A chemical available in the lab.'''
    m_def = Section(label_quantity='formula', categories=[ElnExampleCategory])

    chemical_name = Quantity(type=str, a_eln=dict(component='StringEditQuantity'))

    formula = Quantity(
        type=str,
        description='Empirical formula of the chemical (Hill notation).',
        a_eln=dict(component='StringEditQuantity'))

    form = Quantity(
        type=MEnum([
            'solid',
            'powder',
            'liquid',
            'gas']),
        shape=[],
        description='Physical state of the chemical.',
        a_eln=dict(component='RadioEnumEditQuantity'))

    supplier = Quantity(
        type=str,
        description='Supplier of the chemical.',
        a_eln=dict(component='StringEditQuantity'))

    sku_no = Quantity(
        type=str,
        description=' Stock keeping unit of the chemical. (e.g. 244619-50G)',
        a_eln=dict(component='StringEditQuantity'))

    opening_date = Quantity(
        type=Datetime,
        description='Opening date of the chemical.',
        a_eln=dict(component='DateTimeEditQuantity'))

    impurities = Quantity(
        type=str,
        description='''Descriptions of the impurities of the product.''',
        a_eln=dict(component='StringEditQuantity'))

    cas_number = Quantity(
        type=str,
        description='''The CAS number is a unique and unambiguous
                    identifier for a specific substance.''',
        a_eln=dict(component='StringEditQuantity'))

    sds_link = Quantity(
        type=str,
        description='''The corresponding Safety Data Sheet (SDS) of the product.''',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    comments = Quantity(
        type=str,
        description='''Remarks about the chemical beyond the
                    typically collected information. Here you can collect information
                    about your experience using this chemical, for example if it might
                    have been contaminated during an experiment''',
        a_eln=dict(component='RichTextEditQuantity'))

    def normalize(self, archive, logger):
        super(Chemical, self).normalize(archive, logger)

        if self.formula:
            archive.metadata.entry_name = self.formula
            if self.chemical_name:
                archive.metadata.entry_name += f' ({self.chemical_name})'


class Maintenance(MSection):
    m_def = Section(label_quantity='datetime')

    datetime = Quantity(
        type=Datetime,
        description='The date and time of the maintenance.',
        a_eln=dict(component='DateTimeEditQuantity'))

    maintainer = Quantity(
        type=MEnum([
            'Markus Scheidgen',
            'Pepe Marquez',
            'Sandor Brockhauser',
            'Sherjeel Shabih',
            'Mohammad Nakhaee',
            'David Sitker']),
        description='Name or alias of the persen that performed the maintenance.',
        a_eln=dict(component='AutocompleteEditQuantity'))

    description = Quantity(
        type=str,
        description='Description of what was done to the instrument.',
        a_eln=dict(component='RichTextEditQuantity'))


class Instrument(EntryData):
    '''An instrument available in the lab.'''
    m_def = Section(categories=[ElnExampleCategory])

    name = Quantity(
        type=str,
        description='Name of the instrument used',
        a_eln=dict(component='StringEditQuantity'))

    # Rich text material for ELNs
    description = Quantity(
        type=str,
        description=(
            '''Description of the instrument. May include images, free text
            and tables'''),
        a_eln=dict(component='RichTextEditQuantity'))

    user_manual = Quantity(
        type=str,
        description='''Link to pdf of the user manual''',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    certificate = Quantity(
        type=str,
        description='''Link to pdf of the instrument's certificates''',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    maintenance = SubSection(section_def=Maintenance, repeats=True)

    def normalize(self, archive, logger):
        super(Instrument, self).normalize(archive, logger)

        if self.name:
            archive.metadata.entry_name = self.name


class Process(ArchiveSection):
    ''' Any physical process applied to the sample. '''
    operator = Quantity(
        type=MEnum([
            'Markus Scheidgen',
            'Pepe Marquez',
            'Sandor Brockhauser',
            'Sherjeel Shabih',
            'Mohammad Nakhaee',
            'David Sitker']),
        shape=[],
        description='Name or alias of the process operator.',
        a_eln=dict(
            component='AutocompleteEditQuantity'))

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.',
        a_eln=dict(component='DateTimeEditQuantity'))

    instrument = Quantity(
        type=Reference(Instrument.m_def),
        descriptions='The instrument used for this process.',
        a_eln=dict(component='ReferenceEditQuantity'))

    chemicals = Quantity(
        type=Reference(Chemical.m_def),
        shape=['*'],
        descriptions='The chemicals used in this process',
        a_eln=dict(component='ReferenceEditQuantity'))

    creates_layer = Quantity(
        type=bool,
        a_eln=dict(component='BoolEditQuantity'))

    comments = Quantity(
        type=str,
        description='''Remarks about the process that cannot be seen from the data.
                    Might include rich text, images and potentially tables''',
        a_eln=dict(component='RichTextEditQuantity'))

    def normalize(self, archive, logger):
        super(Process, self).normalize(archive, logger)

        if self.creates_layer:
            logger.debug('create layer if necessary')
            layer_exists = False
            for layer in archive.data.layers:
                if layer.layer_origin is None:
                    continue
                if layer.layer_origin.m_resolved() == self:
                    layer_exists = True

            if not layer_exists:
                logger.debug('create layer')
                archive.data.layers.append(Layer(
                    layer_origin=self,
                    layer_creation_datetime=self.datetime))


class PVDEvaporation(Process):
    '''The physical vapor deposition (PVD) of a layer by evaporation.'''

    def derive_n_values(self):
        if self.process_time is not None:
            return len(self.process_time)
        if self.set_substrate_temperature is not None:
            return len(self.set_substrate_temperature)
        else:
            return 0

    n_values = Quantity(
        type=int, derived=derive_n_values,
        description='Number of registered time values')

    process_time = Quantity(
        type=np.dtype(np.float64), unit='seconds', shape=['n_values'],
        description='The temperature set in the substrate heater')

    set_substrate_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=['n_values'],
        description='The temperature set in the substrate heater',
        a_plot={
            'x': 'process_time', 'y': 'set_substrate_temperature'
        })

    substrate_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=['n_values'],
        description='The temperature measured in the substrate holder',
        a_plot={
            'x': 'process_time', 'y': 'substrate_temperature'
        })

    chamber_pressure = Quantity(
        type=np.dtype(np.float64), unit='pascal', shape=['n_values'],
        description='Data array of the values of the pressure during the process',
        a_plot={
            'x': 'process_time', 'y': 'chamber_pressure'
        })

    number_crucibles = Quantity(
        type=np.dtype(np.int64),
        description='Number of crucibles active in the chamber',
        a_eln=dict(component='NumberEditQuantity', props=dict(maxValue=4, minValue=0)))

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    def normalize(self, archive, logger):
        super(PVDEvaporation, self).normalize(archive, logger)

        if not self.data_file:
            return

        from .pvd import PVDImporter
        importer = PVDImporter()
        with archive.m_context.raw_file(self.data_file) as f:
            self.process_time, substrate_temperature, \
                set_substrate_temperature, chamber_pressure = importer.read(f)
            self.substrate_temperature = substrate_temperature * ureg.degC
            self.set_substrate_temperature = set_substrate_temperature * ureg.degC
            self.chamber_pressure = chamber_pressure * ureg('mbar')


class Targets(MSection):
    '''Targets available in the lab for PLD or Sputtering deposition.'''

    target_name = Quantity(
        type=str,
        description='''Description of the target including information of its
                     chemistry nature . Examples: *Al doped ZnO* or simply *BaZrO3*''')

    elements_in_target = Quantity(type=str, shape=["*"])

    target_diameter = Quantity(type=np.dtype(np.float64), unit='meter')

    supplier = Quantity(
        type=str,
        description='Supplier of the target.')

    opening_date = Quantity(
        type=Datetime,
        description='Opening date of the chemical.')

    impurities = Quantity(
        type=str,
        description='''Descriptions of the impurities of the product.''')

    cas_number = Quantity(
        type=str,
        description='''The CAS number is a unique and unambiguous
                    identifier for a specific substance.''')

    sds_link = Quantity(
        type=str,
        description='''A link to the corresponding Safety Data Sheet (SDS)
                    of the product.''')

    comments = Quantity(
        type=str,
        description='''Remarks about the target beyond the
                    typically collected information in the otehr metadata fields.
                    Here you can collect information about your experience using
                    this target, or additional observations''')


class PLDDeposition(Process):
    '''The deposition process of a layer by Pulsed Laser Deposition (PLD) method.'''

    targets = SubSection(section_def=Targets, repeats=True)


class EbeamEvaporation(Process):
    '''The physical vapor deposition (PVD) of a layer by e-beam evaporation.'''


class HotPlateAnnealing(Process):

    hotplate_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=[],
        description='The temperature set for the hot plate.',
        a_eln=dict(component='NumberEditQuantity'))

    annealing_time = Quantity(
        type=np.dtype(np.float64), unit='seconds', shape=[],
        description='Time of the sample on the hot plate.',
        a_eln=dict(component='NumberEditQuantity'), props=dict(minValue=0))

    relative_humidity = Quantity(
        type=np.dtype(np.float64),  # unit='',
        description=('''
            Relative humidity of the atmosphere in which the experiment was performed.
        '''),
        a_eln=dict(component='NumberEditQuantity', props=dict(minValue=0, maxValue=1)))

    instrument_atmosphere = Quantity(
        type=MEnum([
            'humidity chamber',
            'glove box',
            'ambient']),
        shape=[],
        description='Location or atmosphere in which the process was conducted.',
        a_eln=dict(
            label='Process atmosphere',
            component='RadioEnumEditQuantity'))


class TubeFurnaceAnnealing(Process):
    pass


class RTPAnnealing(Process):
    pass


class SpinCoating(Process):
    pass


class ChemicalBathDeposition(Process):
    pass


class Processes(MSection):
    '''
    Experiment event which  generally change the sample or a new component is added to it.
    For example, in the context of thin films, cleaning the substrate or
    the deposition of a new layer by evaporation are `processes`.
    '''
    m_def = Section(a_eln=dict())

    pvd_evaporation = SubSection(section_def=PVDEvaporation)
    pld_deposition = SubSection(section_def=PLDDeposition)
    ebeam_evaporation = SubSection(section_def=EbeamEvaporation)
    hotplate_annealing = SubSection(section_def=HotPlateAnnealing)
    tubefurnace_annealing = SubSection(section_def=TubeFurnaceAnnealing)
    rtp_annealing = SubSection(section_def=RTPAnnealing)
    spin_coating = SubSection(section_def=SpinCoating)
    chemical_bath_deposition = SubSection(section_def=ChemicalBathDeposition)


class Measurement(ArchiveSection):
    '''
    Any measurement performed on the sample.
    '''
    instrument = Quantity(
        type=Reference(Instrument.m_def),
        descriptions='The instrument used for this measurement.',
        a_eln=dict(component='ReferenceEditQuantity'))

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))


class XrayFluorescence(Measurement):
    '''
    X-ray fluorescence is a technique typically used to obtain information about the
    chemical composition of a sample.
    '''


class XrayDiffraction(Measurement):
    '''
    X-ray diffraction is a technique typically used to characterize the structural
    properties of crystalline materials. The data contains `two_theta` values of the scan
    the corresponding counts collected for each channel
    '''

    start_position_x = Quantity(
        type=np.dtype(np.float64), shape=['*'], unit='meter',
        description='''Length from the left substrate margin to the spot of
                    the first scan. The sample should be oriented with the
                    label at the bottom left corner''',
        a_eln=dict(
            component='NumberEditQuantity',
            label='Length from the left substrate margin'))

    start_position_y = Quantity(
        type=np.dtype(np.float64), shape=['*'], unit='meter',
        description='''Length from the bottom substrate margin to the spot of
                    the first scan. The sample should be oriented with the
                    label at the bottom left corner.''',
        a_eln=dict(
            component='NumberEditQuantity',
            label='Length from the bottom substrate margin'))

    def derive_n_values(self):
        if self.intensity is not None:
            return len(self.intensity)
        if self.two_theta is not None:
            return len(self.two_theta)
        else:
            return 0

    n_values = Quantity(type=int, derived=derive_n_values)

    two_theta = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        unit='radian',
        description='The 2-theta range of the difractogram',
        a_plot={
            'x': 'two_theta', 'y': 'intensity'
        })

    q_vector = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        unit='meter**(-1)',
        description='The scattering vector *Q* of the difractogram',
        a_plot={
            'x': 'q_vector', 'y': 'intensity'
        })

    intensity = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        description='The count at each 2-theta value, dimensionless',
        a_plot={
            'x': 'two_theta', 'y': 'intensity'
        })

    omega = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        unit='radian',
        description='The omega range of the difractogram')

    xray_tube_current = Quantity(
        type=np.dtype(np.float64),
        unit='A',
        description='Current of the X-ray tube',
        a_eln=dict(
            component='NumberEditQuantity',
            label='Current of the X-ray tube',
            props=dict(minValue=1e4, maxValue=4e4)))

    xray_tube_voltage = Quantity(
        type=np.dtype(np.float64),
        unit='V',
        description='Voltage of the X-ray tube',
        a_eln=dict(
            component='NumberEditQuantity',
            label='Voltage of the X-ray tube',
            props=dict(minValue=1e4, maxValue=4e4)))

    kalpha_one = Quantity(
        type=np.dtype(np.float64),
        unit='angstrom',
        description='Wavelength of the Kα1 line')

    kalpha_two = Quantity(
        type=np.dtype(np.float64),
        unit='angstrom',
        description='Wavelength of the Kα2 line')

    ratio_kalphatwo_kalphaone = Quantity(
        type=np.dtype(np.float64),
        description='Kα2/Kα1 intensity ratio')

    kbeta = Quantity(
        type=np.dtype(np.float64),
        unit='angstrom',
        description='Wavelength of the Kβ line')

    scan_axis = Quantity(
        type=str,
        description='Axis scanned')

    integration_time = Quantity(
        type=np.dtype(np.float64),
        shape=['*'],
        description='Integration time per channel')

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'))

    def normalize(self, archive, logger):
        super(XrayDiffraction, self).normalize(archive, logger)

        if not self.data_file:
            return

        import xrdtools

        with archive.m_context.raw_file(self.data_file) as f:
            xrdml_dict = xrdtools.read_xrdml(f.name)
            self.intensity = xrdml_dict['data']
            self.two_theta = xrdml_dict['x'] * ureg('degree')
            self.omega = xrdml_dict['Omega'] * ureg('degree')
            self.kalpha_one = xrdml_dict['kAlpha1']
            self.kalpha_two = xrdml_dict['kAlpha2']
            self.ratio_kalphatwo_kalphaone = xrdml_dict['kAlphaRatio']
            self.kbeta = xrdml_dict['kBeta']
            self.scan_axis = xrdml_dict['scanAxis']
            self.integration_time = xrdml_dict['time']
            self.q_vector = (4 * np.pi / self.kalpha_one) * np.sin((self.two_theta))


class RamanSpectroscopy(Measurement): pass


class Resistivity(Measurement): pass


class SolarSimulator(Measurement): pass


class CapacitanceVoltage(Measurement): pass


class EQE(Measurement): pass


class SteadyStatePL(Measurement): pass


class PLImaging(Measurement): pass


class TimeResolvedPL(Measurement): pass


class UVVisNIRImaging(Measurement): pass


class UVVisNIRSpectroscopy(Measurement): pass


class TerahertzSpectroscopy(Measurement): pass


class Measurements(MSection):
    '''
    Experimental procedure in which a sample gets characterized
    by a technique. For example, a measurment by X-ray diffraction to characterize
    the structural properties of an specimen or X-ray fluorescence
    to characterize its composition.
    '''
    m_def = Section(a_eln=dict())

    xray_diffraction = SubSection(section_def=XrayDiffraction, repeats=True)
    xray_fluorescence = SubSection(section_def=XrayFluorescence, repeats=True)
    raman_spectroscopy = SubSection(section_def=RamanSpectroscopy, repeats=True)
    resistivity = SubSection(section_def=Resistivity, repeats=True)
    solar_simulator = SubSection(section_def=SolarSimulator, repeats=True)
    capacitance_voltage = SubSection(section_def=CapacitanceVoltage, repeats=True)
    eqe = SubSection(section_def=EQE, repeats=True)
    steady_state_pl = SubSection(section_def=SteadyStatePL, repeats=True)
    pl_imaging = SubSection(section_def=PLImaging, repeats=True)
    time_resolved_pl = SubSection(section_def=TimeResolvedPL, repeats=True)
    uv_vis_nir_imaging = SubSection(section_def=UVVisNIRImaging, repeats=True)
    uv_vis_spectroscopy = SubSection(section_def=UVVisNIRSpectroscopy, repeats=True)
    terahertz_spectroscopy = SubSection(section_def=TerahertzSpectroscopy, repeats=True)


class DerivedData(MSection):
    '''
    Additional information that is gained from already existing data i.e. from `processes`,
     `measurements`  and/or other `derived_data`. Examples for these are smoothing of data,
    background subtraction, calculations (e.g. transmittance from raw spectra,
    bandgap from the absorption coefficient) or simulations.
    '''


class PhysicalProperties(MSection):
    '''
    It describes physical properties of the layer, like its dimensions in the *x* and *y*
    directions or its thickness.
    '''
    m_def = Section(a_eln=dict())

    x_length = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=[],
        description='Dimension of the layer in the *x*-direction',
        a_eln=dict(component='NumberEditQuantity'))

    y_length = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=[],
        description='Dimension of the layer in the *y*-direction',
        a_eln=dict(component='NumberEditQuantity'))

    thickness = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=["*"],
        description='Thickness of the layer',
        a_eln=dict(component='NumberEditQuantity'))


class CompositionalProperties(MSection):
    '''
    It describes the chemical properties of the layer, like the *elements* contained
    in the layer or the *compounds* if they are known. These properties can be filled with
    knowledge obtained from `measurements`.
    '''
    elements = Quantity(type=str, shape=["*"])
    chemical_formula = Quantity(type=str)


class StructuralProperties(MSection): pass


class OptoelectronicProperties(MSection): pass


class Layer(MSection):
    m_def = Section(a_eln=dict())

    layer_type = Quantity(
        type=MEnum([
            'substrate',
            'absorber layer',
            'metal contact',
            'buffer layer',
            'p-type contact',
            'n-type contact']),
        shape=[],
        description='type of the layer',
        a_eln=dict(component='EnumEditQuantity'))

    layer_creation_datetime = Quantity(
        type=Datetime,
        description='Creation date of the layer.',
        a_eln=dict(
            label='Creation date and time of the layer.',
            component='DateTimeEditQuantity'))

    layer_position = Quantity(
        type=np.dtype(np.int64),
        description='''Position of the layer within
        the stack counting from substrate = 0; not necessarily gapless''',
        a_eln=dict(component='NumberEditQuantity'))

    layer_origin = Quantity(
        type=Reference(Process.m_def),
        description='A link to the `process` where the layer was created.')

    x_length = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=[],
        description='Dimension of the layer in the *x*-direction',
        a_eln=dict(component='NumberEditQuantity'))

    y_length = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=[],
        description='Dimension of the layer in the *y*-direction',
        a_eln=dict(component='NumberEditQuantity'))

    thickness = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=["*"],
        description='Thickness of the layer',
        a_eln=dict(component='NumberEditQuantity'))

    elements = Quantity(type=str, shape=["*"], a_eln=dict(component='StringEditQuantity'))
    chemical_formula = Quantity(type=str, a_eln=dict(component='StringEditQuantity'))

    # physical_properties = SubSection(section_def=PhysicalProperties)
    # compositional_properties = SubSection(section_def=CompositionalProperties)
    # structural_properties = SubSection(section_def=StructuralProperties)
    # optoelectronic_properties = SubSection(section_def=OptoelectronicProperties)


class Projects(MSection):
    project_acronym = Quantity(
        type=str, description='Acronim or short name of the project.')

    project_long_name = Quantity(
        type=str, description='Long name of the prioject.')

    # Might cotain images and free text in ELNs features
    project_description = Quantity(
        type=str, description='A description of the project activities.')

    project_website = Quantity(
        type=str, description='Project website.')

    funding_agency = Quantity(
        type=str, description='The physical state of the sample.')

    gran_agreement_id = Quantity(
        type=str, description='ID of the contract of that is funding this project.')


class Sample(EntryData):
    '''
    Sample in which one or various properties
    can vary across the area of the substrate,
    thus containing subsamples whithin.
    An example would be a subtrate with a deposited film on top in which
    the chemical compostion of the film varies in the x and y directions.
    '''
    m_def = Section(categories=[ElnExampleCategory])

    sample_owner = Quantity(
        type=MEnum([
            'Markus Scheidgen',
            'Pepe Marquez',
            'Sandor Brockhauser',
            'Sherjeel Shabih',
            'Mohammad Nakhaee',
            'David Sitker']),
        shape=[],
        description='Name or alias of the process operator.',
        a_eln=dict(component='AutocompleteEditQuantity'))

    sample_id = Quantity(
        type=str,
        description='Full sample id.',
        a_eln=dict(component='StringEditQuantity'))

    sample_name = Quantity(
        type=str,
        description='Short name of the sample (the id on the substrate), e.g. `4001-8`.',
        a_eln=dict(component='StringEditQuantity'))

    creation_datetime = Quantity(
        type=Datetime,
        description='Creation date of the sample.',
        a_eln=dict(component='DateTimeEditQuantity'))

    institute = Quantity(
        type=str,
        description='Alias/short name of the home institute of the owner, i.e. `HZB`.',
        a_eln=dict(component='StringEditQuantity'))

    description = Quantity(
        type=str, a_eln=dict(component='RichTextEditQuantity'))

    processes = SubSection(section_def=Processes)
    measurements = SubSection(section_def=Measurements)
    # derived_data = SubSection(section_def=DerivedData, repeats=True)
    layers = SubSection(section_def=Layer, repeats=True)
    # projects = SubSection(section_def=Projects, repeats=True)


m_package.__init_metainfo__()
