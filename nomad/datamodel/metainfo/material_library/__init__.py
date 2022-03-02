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

from matplotlib.afm import CharMetrics
import numpy as np

from nomad.units import ureg
from nomad.metainfo import (
    MSection, Package, Quantity, SubSection, MEnum, Reference, Datetime, Section)
from nomad.metainfo.metainfo import SectionProxy
from nomad.datamodel.data import EntryData


m_package = Package(name='material_library')


class Chemicals(MSection):
    '''Chemicals available in the lab.'''

    formula = Quantity(
        type=str,
        description='Empirical formula of the chemical (Hill notation).')

    form = Quantity(
        type=MEnum([
            'solid',
            'powder',
            'liquid',
            'gas']),
        shape=[],
        description='Physical state of the chemical.')

    supplier = Quantity(
        type=str,
        description='Supplier of the chemical.')

    opening_date = Quantity(
        type=Datetime,
        description='Opening date of the chemical.')

    # pack_size = Quantity(
    #     type=np.float64, unit = ['mg', 'ml'],
    #     description='Distributor of the chemical.')

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
        description='''Remarks about the chemical beyond the
                    typically collected information. Here you can collect information
                    about your experience using this chemical, for exmaple if it might
                    have been contaminated during an experiment''',
        a_eln=dict(component='RichTextEditQuantity'))


class Instrument(MSection):
    '''Chemicals available in the lab.'''

    instrument_name = Quantity(
        type=str,
        description='Name of the instrument used')

    # Rich text material for ELNs
    instrument_information = Quantity(
        type=str,
        description='''Description of the instrument. May include images, free text
                    and tables''',
        a_eln=dict(component='RichTextEditQuantity'))

    user_manual = Quantity(
        type=str,
        description='''Link to pdf of the user manual''')


class PVDEvaporation(MSection):
    '''The physical vapor deposition (PVD) of a layer by evaporation.'''

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
            label='Layer type',
            component='AutocompleteEditQuantity'))

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    comments = Quantity(
        type=str,  # Revise type
        description='''Remarks about the process that cannot be seen from the data.
                    Might include rich text, images and potentially tables''',
        a_eln=dict(component='RichTextEditQuantity'))

    def derive_n_values(self):
        if self.process_time is not None:
            return len(self.process_time)
        if self.set_substrate_temperature is not None:
            return len(self.set_substrate_temperature)
        else:
            return 0

    # Ask Markus about data type and tell that should have the shape of the time array
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
        description='Data array of the values of the pressure during the process')

    number_crucibles = Quantity(
        type=np.dtype(np.int64),  # unit='',
        description='Number of crucibles active in the chamber',
        a_eln=dict(component='NumberEditQuantity', props=dict(maxValue=4, minValue=0)))

    comments = Quantity(
        type=str,  # Revise type
        description='''Remarks about the process that cannot be seen from the data.
                    Might include rich text, images and potentially tables''',
        a_eln=dict(component='RichTextEditQuantity'))

    # layer_created = SubSection(section_def=Layers, repeats=True)
    instrument = SubSection(section_def=Instrument)
    chemicals = SubSection(section_def=Chemicals, repeats=True)

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'))

    def normalize(self, archive, logger):
        from nomad.datamodel.metainfo.material_library.PvdPImporter import Importer
        importer = Importer()
        if (self.data_file):
            with archive.m_context.raw_file(self.data_file) as f:
                self.process_time, substrate_temperature, \
                    set_substrate_temperature, chamber_pressure = importer.read(f)
                self.substrate_temperature = substrate_temperature * ureg.degC
                self.set_substrate_temperature = set_substrate_temperature * ureg.degC
                self.chamber_pressure = chamber_pressure * ureg('mbar')
                # self.n_values = self.n_values


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

    # pack_size = Quantity(
    #     type=np.float64, unit = ['mg', 'ml'],
    #     description='Distributor of the chemical.')

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


class PLDDeposition(MSection):
    '''The deposition process of a layer by Pulsed Laser Deposition (PLD) method.'''

    operator = Quantity(
        type=str,
        description='Name or alias of the process operator.')

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    targets = SubSection(section_def=Targets, repeats=True)

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.')

    instrument = SubSection(section_def=Instrument)
# class ProcessEntries(MSection):
#     '''The physical vapor deposition (PVD) of a layer by e-beam evaporation.'''

#     operator = Quantity(
#         type=str,
#         description='Name or alias of the process operator.')

#     datetime = Quantity(
#         type=Datetime,
#         description='Finishing date and time of the process.')

#     comments = Quantity(
#         type=str,
#         description='Remarks about the process that cannot be seen from the data.')


class EbeamEvaporation(MSection):
    '''The physical vapor deposition (PVD) of a layer by e-beam evaporation.'''

    operator = Quantity(
        type=str,
        description='Name or alias of the process operator.')

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.')


class HotPlateAnnealing(MSection):

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
            label='Layer type',
            component='AutocompleteEditQuantity'))

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    hotplate_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=[],
        description='The temperature set in the hot plate',
        a_eln=dict(component='NumberEditQuantity'))

    annealing_time = Quantity(
        type=np.dtype(np.float64), unit='seconds', shape=[],
        description='Time of the sample on the hot plate',
        a_eln=dict(component='NumberEditQuantity', props=dict(minValue=0, maxValue=1)))

    relative_humidity = Quantity(
        type=np.dtype(np.float64),  # unit='',
        description='''Relative humidity of the atmosphere in which the experiment was
                    performed''',
        a_eln=dict(component='NumberEditQuantity', props=dict(minValue=0, maxValue=1)))

    atmosphere = Quantity(
        type=MEnum([
            'humidity chamber',
            'glove box',
            'N2',
            'Ar',
            'Ambient']),
        shape=[],
        description='Atmosphere in which the process was conducted.',
        a_eln=dict(
            label='Process atmosphere',
            component='AutocompleteEditQuantity'))

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))


class TubeFurnaceAnnealing(MSection):

    operator = Quantity(
        type=str,
        description='Name or alias of the process operator.')

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))


class RTPAnnealing(MSection):

    operator = Quantity(
        type=str,
        description='Name or alias of the process operator.')

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))


class SpinCoating(MSection):

    operator = Quantity(
        type=str,
        description='Name or alias of the process operator.')

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))


class ChemicalBathDeposition(MSection):

    operator = Quantity(
        type=str,
        description='Name or alias of the process operator.')

    datetime = Quantity(
        type=Datetime,
        description='Finishing date and time of the process.')

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))


class Processes(MSection):
    '''
    Experimetal event which  generally change the library or a new component is added to it.
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


class XrayFluorescence(MSection):
    '''
    X-ray fluorescence is a technique typically used to obtain information about the
    chemical composition of a sample.
    '''


# TODO include mapping of positions for several imported files
class XrayDiffraction(MSection):
    '''
    X-ray diffraction is a technique typically used to characterize the structural
    properties of crystalline materials. The data contains `two_theta` values of the scan
    the corresponding counts collected for each channel
    '''

    # def derive_x_values(self):
    #     if self.count is not None:
    #         return len(self.x_positions)
    #     else:
    #         return 0

    # x_values = Quantity(type=int, derived=derive_x_values)

    # def derive_y_values(self):
    #     if self.count is not None:
    #         return len(self.y_positions)
    #     else:
    #         return 0

    # y_values = Quantity(type=int, derived=derive_y_values)

    # start_position_x = Quantity(
    #     type=np.dtype(np.float64), shape=['*'], unit='meter',
    #     description='''Length from the left substrate margin to the spot of
    #                 the first scan. The sample should be oriented with the
    #                 label at the bottom left corner''')

    # start_position_y = Quantity(
    #     type=np.dtype(np.float64), shape=['*'], unit='meter',
    #     description='''Length from the bottom substrate margin to the spot of
    #                 the first scan. The sample should be oriented with the
    #                 label at the bottom left corner.''')

    # x_positions = Quantity(
    #     type=np.dtype(np.float64), shape=['x_values'], unit='meter',
    #     description='''Scan positions on the *x* direction within the substrates.
    #                 The origin of coordinates is at the labelled corner of the substrate
    #                 from the top view ''')

    # y_positions = Quantity(
    #     type=np.dtype(np.float64), shape=['y_values'], unit='meter',
    #     description='''Scan positions on the *y* direction within the substrates
    #                 The origin of coordinates is at the labelled corner of the substrate
    #                 from the top view ''')

    def derive_n_values(self):
        if self.intensity is not None:
            return len(self.intensity)
        if self.two_theta is not None:
            return len(self.two_theta)
        else:
            return 0

    n_values = Quantity(type=int, derived=derive_n_values)

    # two_theta = Quantity(
    #     type=np.dtype(np.float64), shape=['x_values', 'y_values', 'n_values'],
    #     unit='degrees',
    #     description='The 2-theta range of the difractogram')

    # intensity = Quantity(
    #     type=np.dtype(np.float64), shape=['x_values', 'y_values', 'n_values'],
    #     description='The count at each 2-theta value, dimensionless',
    #     a_plot={
    #         'x': 'two_theta', 'y': 'intensity'
    #     })

    two_theta = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        unit='radian',
        description='The 2-theta range of the difractogram')

    intensity = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        description='The count at each 2-theta value, dimensionless',
        a_plot={
            'x': 'two_theta', 'y': 'intensity'
        })

    xray_tube_material = Quantity(
        type=str,
        shape=["*"],
        description='Material of the X-ray tube')

    xray_tube_current = Quantity(
        type=np.dtype(np.float64),
        unit='A',
        description='Current of the X-ray tube')

    xray_tube_voltage = Quantity(
        type=np.dtype(np.float64),
        unit='V',
        description='Voltage of the X-ray tube')

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
        shape=['*'],
        description='Axis scanned')

    integration_time = Quantity(
        type=np.dtype(np.float64),
        description='Integration time per channel')

    instrument = SubSection(section_def=Instrument)

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'))

    comments = Quantity(
        type=str,
        description='Remarks about the process that cannot be seen from the data.',
        a_eln=dict(component='RichTextEditQuantity'))

    # def normalize(self, archive, logger):
    #     from nomad.datamodel.metainfo.material_library.XRDImporter import Importer
    #     importer = Importer()
    #     if (self.data_file):
    #         with archive.m_context.raw_file(self.data_file) as f:
    #             self.intensity,
    #             self.two_theta,
    #             self.kalpha_one,
    #             self.kalpha_two,
    #             self.ratio_kalphatwo_kalphaone,
    #             self.kbeta, self.scan_axis,
    #             self.integration_time = importer.read(f)

    def normalize(self, archive, logger):
        import xrdtools
        if (self.data_file):
            with archive.m_context.raw_file(self.data_file) as f:
                xrdml_dict = xrdtools.read_xrdml(f.name)
                self.intensity = xrdml_dict['data']
                self.two_theta = xrdml_dict['x'] * ureg('°')
                self.kalpha_one = xrdml_dict['kAlpha1']
                self.kalpha_two = xrdml_dict['kAlpha2']
                self.ratio_kalphatwo_kalphaone = xrdml_dict['kAlphaRatio']
                self.kbeta = xrdml_dict['kBeta']
                self.scan_axis = xrdml_dict['scanAxis']
                self.integration_time = xrdml_dict['time']


class RamanSpectroscopy(MSection):
    '''

    '''


class Resistivity(MSection):
    '''

    '''


class SolarSimulator(MSection):
    '''

    '''


class CapacitanceVoltage(MSection):
    '''

    '''


class EQE(MSection):
    '''

    '''


class SteadyStatePL(MSection):
    '''

    '''


class PLImaging(MSection):
    '''

    '''


class TimeResolvedPL(MSection):
    '''

    '''


class UVVisNIRImaging(MSection):
    '''

    '''


class UVVisNIRSpectroscopy(MSection):
    '''

    '''


class TerahertzSpectroscopy(MSection):
    '''

    '''


class Measurements(MSection):
    '''
    Experimental procedure in which a material library or sample gets characterized
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
        description='Dimension of the layer in the *x*-direction')

    y_length = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=[],
        description='Dimension of the layer in the *y*-direction')

    thickness = Quantity(
        type=np.dtype(np.float64), unit='meter', shape=["*"],
        description='Thickness of the layer')


class CompositionalProperties(MSection):
    '''
    It describes the chemical properties of the layer, like the *elements* contained
    in the layer or the *compounds* if they are known. These properties can be filled with
    knowledge obtained from `measurements`.
    '''
    elements = Quantity(type=str, shape=["*"])
    chemical_formula = Quantity(type=str)


class StructuralProperties(MSection):
    '''

    '''


class OptoelectronicProperties(MSection):
    '''

    '''


class Layers(MSection):
    '''
    List of layers contained in the library, which in turn have subsections
    describing the layer and its properties.
    '''
    m_def = Section(a_eln=dict())

    layer_id = Quantity(
        type=np.dtype(np.int64), description='''unique identifier of the layer
        (counting in the order of adding).
        "0" would be assigned to the substrate.''')

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
        a_eln=dict(
            label='Layer type',
            component='EnumEditQuantity'))

    layer_creation_datetime = Quantity(
        type=Datetime,
        description='Creation date of the layer.')

    layer_position = Quantity(
        type=np.dtype(np.int64),
        description='''Position of the layer within
        the stack counting from substrate = 0; not necessarily gapless''')

    layer_origin = Quantity(
        type=str, description='A link to the `process` where the layer was created.')

    physical_properties = SubSection(section_def=PhysicalProperties)
    compositional_properties = SubSection(section_def=CompositionalProperties)
    structural_properties = SubSection(section_def=StructuralProperties)
    optoelectronic_properties = SubSection(section_def=OptoelectronicProperties)


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


class MaterialLibrary(EntryData):
    '''
    Sample in which one or various properties
    can vary across the area of the substrate,
    thus containing subsamples whithin.
    An example would be a subtrate with a deposited film on top in which
    the chemical compostion of the film varies in the x and y directions.
    '''

    library_id = Quantity(
        type=str,
        description='Full library id.')

    library_name = Quantity(
        type=str,
        description='Short name of the library (the id on the substrate), e.g. `4001-8`.',
        a_eln=dict(component='StringEditQuantity'))

    creation_datetime = Quantity(
        type=Datetime,
        description='Creation date of the library.',
        a_eln=dict(component='DateTimeEditQuantity'))

    institute = Quantity(
        type=str,
        description='Alias/short name of the home institute of the owner, i.e. `HZB`.',
        a_eln=dict(component='StringEditQuantity'))

    description = Quantity(
        type=str, a_eln=dict(component='RichTextEditQuantity'))

    processes = SubSection(section_def=Processes)
    measurements = SubSection(section_def=Measurements)
    derived_data = SubSection(section_def=DerivedData, repeats=True)
    layers = SubSection(section_def=Layers, repeats=True)
    projects = SubSection(section_def=Projects, repeats=True)


m_package.__init_metainfo__()
