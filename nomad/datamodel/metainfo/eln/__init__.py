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
from nomad.metainfo import Package, Quantity, Datetime, Reference, Section, SubSection
from ase.data import chemical_symbols, atomic_numbers, atomic_masses
from nomad.datamodel.metainfo.eln.perovskite_solar_cell_database import (
    add_solar_cell, add_band_gap
)

from nomad.datamodel.metainfo.eln.nexus_data_converter import (
    NexusDataConverter,
    ElnYamlConverter
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
            data. Will be generated automatically upon saving if name, users and institute
            is filled.''',
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


class Entity(ElnBaseSection):
    ''' A ELN base section that can be used for chemicals.'''
    pass


class Activity(ElnBaseSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used for
    laboratory activities, e.g. processes, characterizations, measurements, etc.
    '''
    datetime = Quantity(
        type=Datetime,
        description='The date and time when this activity was started.',
        a_eln=dict(component='DateTimeEditQuantity', label='Starting Time'))

    end_time = Quantity(
        type=Datetime,
        description='The date and time when this activity was done.',
        a_eln=dict(component='DateTimeEditQuantity', label='Ending Time'))

    method = Quantity(
        type=str,
        description='A short consistent handle for the applied method.')

    location = Quantity(
        type=str,
        description='location of the experiment.',
        a_eln=dict(component='StringEditQuantity'))

    def normalize(self, archive, logger):
        super(Activity, self).normalize(archive, logger)

        if archive.results.eln.methods is None:
            archive.results.eln.methods = []
        if self.method:
            archive.results.eln.methods.append(self.method)
        else:
            archive.results.eln.methods.append(self.m_def.name)


class ElementalComposition(MSection):
    '''A section for describing the elemental composition of a system, i.e. the element
    and its atomic fraction.
    '''
    m_def = Section(label_quantity='element')
    element = Quantity(
        type=MEnum(chemical_symbols[1:]),
        description='''
        The symbol of the element, e.g. 'Pb'.
        ''',
        a_eln=dict(component='EnumEditQuantity'))
    atomic_fraction = Quantity(
        type=np.float64,
        description='''
        The atomic fraction of the element in the system it is contained within.
        Per definition a positive value less than or equal to 1.
        ''',
        a_eln=dict(component='NumberEditQuantity')
    )


class System(Entity):
    '''A base class for any system of materials which is investigated or used to construct
    other systems.'''
    elemental_composition = SubSection(
        description='''
        A list of all the elements found in the system together and their respective
        atomic fraction within the system.
        ''',
        section_def=ElementalComposition,
        repeats=True)

    def normalize(self, archive, logger: Any) -> None:
        '''The normalizer for the `System` class. Will add a results.material subsection
        if none exists. Will populate the elements property of that subsection with a
        uniques list of element symbols retrieved from the system elemental_composition.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (Any): A structlog logger.
        '''
        super(System, self).normalize(archive, logger)
        if not archive.results.material:
            archive.results.material = Material()
        elements = set()
        for comp in self.elemental_composition:
            if comp.element not in chemical_symbols:
                logger.warn(
                    f"'{comp.element}' is not a valid element symbol and this "
                    "elemental_composition section will be ignored.")
            else:
                elements.add(comp.element)
        archive.results.material.elements = list(elements)


class Instrument(ElnBaseSection):
    ''' A ELN base section that can be used for instruments.'''
    def normalize(self, archive, logger):
        super(Instrument, self).normalize(archive, logger)

        if self.name:
            if archive.results.eln.instruments is None:
                archive.results.eln.instruments = []
            archive.results.eln.instruments.append(self.name)


class Process(Activity):
    ''' A ELN base section that can be used for processes.'''
    pass


class Measurement(Activity):
    ''' A ELN base section that can be used for measurements.'''
    pass


class Component(ArchiveSection):
    '''A section for describing a component and its role in an ensemble.'''
    name = Quantity(
        type=str,
        description='A short name for the component.',
        a_eln=dict(component='StringEditQuantity', label='Component name'))
    system = Quantity(
        type=Reference(System.m_def),
        description='A reference to the component system.',
        a_eln=dict(component='ReferenceEditQuantity'))
    mass = Quantity(
        type=np.float64,
        description='The mass of the component.',
        unit='kg',
        a_eln=dict(component='NumberEditQuantity', defaultDisplayUnit='mg'))

    def normalize(self, archive, logger: Any) -> None:
        '''The normalizer for the `Component` class. If none is set, the normalizer
        will set the name of the component to be that of the referenced system if it has
        one.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (Any): A structlog logger.
        '''
        super(Component, self).normalize(archive, logger)
        if self.name is None and self.system is not None:
            self.name = self.system.name


class Ensemble(System):
    '''A base class for an ensemble of material systems. Each component of the ensemble is
    of a (sub)type of `System`.'''
    components = SubSection(
        description='''
        A list of all the components of the ensemble containing a name, reference to the
        system section and mass of that component.
        ''',
        section_def=Component,
        repeats=True)

    @staticmethod
    def _atomic_to_mass(composition: List[ElementalComposition], mass: float) -> Dict[str, float]:
        '''Private static method for converting list of ElementalComposition objects to
        dictionary of element masses with the element symbol as key and mass as value.

        Args:
            composition (List[ElementalComposition]): List of ElementalComposition objects
            containing the element symbol and the atomic fraction of that element
            respectively.
            mass (float): The total mass of that component which is simply multiplied with
            the atomic fractions.

        Returns:
            Dict[str, float]: Dictionary of element masses with the element symbol as key
            and the mass as (float) value. If any atomic fraction is None, all masses will
            be None.
        '''
        atomic_fractions = []
        atom_masses = []
        for comp in composition:
            if None in [mass, comp.element, comp.atomic_fraction]:
                return {comp.element: None for comp in composition if comp.element}
            atomic_fractions.append(comp.atomic_fraction)
            if comp.element not in atomic_numbers:
                raise ValueError(f'Unknown element symbol: {comp.element}')
            atom_masses.append(atomic_masses[atomic_numbers[comp.element]])
        masses = np.array(atomic_fractions) * np.array(atom_masses)
        mass_array = mass * masses / masses.sum()
        mass_dict: Dict[str, float] = {}
        for c, m in zip(composition, mass_array):
            if c.element in mass_dict:
                mass_dict[c.element] += m
            else:
                mass_dict[c.element] = m
        return mass_dict

    @staticmethod
    def _mass_to_atomic(mass_dict: Dict[str, float]) -> List[ElementalComposition]:
        '''Private static method for converting ditctionary of elements with their masses
        to a list of ElementalComposition objects containing atomic fractions.

        Args:
            mass_dict (Dict[str, float]): Dictionary of element masses with the element
            symbol as key and the mass as (float) value.

        Returns:
            List[ElementalComposition]: List of ElementalComposition objects containing
            the element symbol and the atomic fraction of that element respectively. If
            any input mass is None all the "atomic_fraction" properties are set to None.
        '''
        atoms = []
        elements = []
        for element, mass in mass_dict.items():
            elements.append(element)
            if mass is None:
                return [ElementalComposition(element=symbol) for symbol in mass_dict]
            else:
                atoms.append(mass / atomic_masses[atomic_numbers[element]])
        atoms_array: np.ndarray = np.array(atoms).astype(float)
        total_atoms = atoms_array.sum()
        atomic_fractions = atoms_array / total_atoms

        return [ElementalComposition(element=symbol, atomic_fraction=atomic_fraction)
                for atomic_fraction, symbol in zip(atomic_fractions, elements)]

    def normalize(self, archive, logger: Any) -> None:
        '''The normalizer for the `Ensemble` class. If the elemental composition list is
        empty, the normalizer will iterate over the components and extract all the
        elements for populating the elemental composition list. If masses are provided for
        all components and the elemental composition of all components contain atomic
        fractions the normalizer will also calculate the atomic fractions for the
        ensemble. The populated elemental composition list is added to the results by the
        normalizer in the `System` super class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (Any): A structlog logger.
        '''
        if logger is None:
            logger = utils.get_logger(__name__)
        if not self.elemental_composition:
            from pint import UnitRegistry
            ureg = UnitRegistry()

            mass_dict: Dict[str, float] = {}
            for component in self.components:
                if component.mass is None:
                    mass = None
                else:
                    mass = component.mass.to(ureg.kilogram).magnitude
                component_dict = self._atomic_to_mass(
                    component.system.elemental_composition, mass)
                for element, mass in component_dict.items():
                    if element in mass_dict and mass is not None:
                        mass_dict[element] += mass
                    else:
                        mass_dict[element] = mass
            self.elemental_composition = self._mass_to_atomic(mass_dict)

        super(Ensemble, self).normalize(archive, logger)


class CASExperimentalProperty(MSection):
    '''A section for experimental properties retrieved from the CAS API.'''

    name = Quantity(
        type=str,
        description='CAS experimental property name.')

    property = Quantity(
        type=str,
        description='CAS experimental property.')

    sourceNumber = Quantity(
        type=str,
        description='CAS experimental property source.')


class CASPropertyCitation(MSection):
    '''A section for citations of the experimental properties retrieved from the CAS API.
    '''

    docUri = Quantity(
        type=str,
        description='CAS property citation document uri.')

    sourceNumber = Quantity(
        type=int,
        decription='CAS property citation source number.')

    source = Quantity(
        type=str,
        description='CAS property citation source.')


class Substance(System):
    '''A base section for any substance defined in the ELN.'''

    name = Quantity(
        type=str,
        description='The name of the substance entry.',
        a_eln=dict(component='StringEditQuantity', label='Substance name'))

    lab_id = Quantity(
        type=str,
        description='''
        A human human readable substance ID that is at least unique for the lab.
        ''',
        a_eln=dict(component='StringEditQuantity', label='Substance ID'))

    cas_uri = Quantity(
        type=str,
        description='CAS uri',
        a_eln=dict(component='StringEditQuantity', label='CAS uri'))

    cas_number = Quantity(
        type=str,
        description='CAS number.',
        a_eln=dict(component='StringEditQuantity', label='CAS number'))

    cas_name = Quantity(
        type=str,
        description='CAS name.',
        a_eln=dict(component='StringEditQuantity', label='CAS name'))

    image = Quantity(
        type=str,
        description='CAS image.',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor', label='Image of substance'))

    inchi = Quantity(
        type=str,
        description='CAS inchi.',
        a_eln=dict(component='StringEditQuantity'))

    inchi_key = Quantity(
        type=str,
        description='CAS inchi key.',
        a_eln=dict(component='StringEditQuantity'))

    smile = Quantity(
        type=str,
        description='CAS smile.',
        a_eln=dict(component='StringEditQuantity'))

    canonical_smile = Quantity(
        type=str,
        description='CAS canonical smile.',
        a_eln=dict(component='StringEditQuantity'))

    molecular_formula = Quantity(
        type=str,
        description='CAS molecular formula.',
        a_eln=dict(component='StringEditQuantity'))

    molecular_mass = Quantity(
        type=np.dtype(np.float64),
        unit='Da',
        description='CAS molecular mass.',
        a_eln=dict(
            component='NumberEditQuantity'))

    cas_experimental_properties = SubSection(
        section_def=CASExperimentalProperty,
        repeats=True)

    cas_property_citations = SubSection(
        section_def=CASPropertyCitation,
        repeats=True)

    cas_synonyms = Quantity(
        type=str,
        shape=['*'],
        description='CAS synonyms.'
    )

    description = Quantity(
        type=str,
        description='''
        A field for adding additional information about the substance that is not captured
        by the other quantities and subsections.
        ''',
        a_eln=dict(
            component='RichTextEditQuantity', label='Detailed substance description'))

    def _populate_from_cas(self, archive, logger: Any) -> None:
        '''Private method for populating the attributes from a call to the CAS API using
        the `cas_number`.
        Will overwrite exisiting CAS attributes if the query provides a value for them.
        I.e. all attributes that begin with `cas_`.

        Args:
            archive (EntryArchive): The archive that is being normalized.
            logger (Any): A structlog logger.
        '''
        import httpx
        response = httpx.get(
            f'https://commonchemistry.cas.org/api/detail?cas_rn={self.cas_number}',
            timeout=2)
        if response.status_code == 200:
            response_dict = response.json()

            self.cas_uri = response_dict.get("uri", None)

            cas_name = response_dict.get("name", None)
            if cas_name:
                self.cas_name = re.sub(r'<.*?>', '', cas_name)
                if not self.name:
                    self.name = self.cas_name

            if not self.image:
                image = response_dict.get("image", None)
                if image:
                    self.image = f"cas_{self.cas_number}_image.svg"
                    with archive.m_context.raw_file(self.image, 'w') as fh:
                        fh.write(image)

            if not self.inchi:
                self.inchi = response_dict.get("inchi", None)

            if not self.inchi_key:
                self.inchi_key = response_dict.get("inchiKey", None)

            if not self.smile:
                self.smile = response_dict.get("smile", None)

            if not self.canonical_smile:
                self.canonical_smile = response_dict.get("canonicalSmile", None)

            if not self.molecular_formula:
                molecular_formula = response_dict.get("molecularFormula", None)
                if molecular_formula:
                    self.molecular_formula = re.sub(r'<.*?>', '', molecular_formula)

            if not self.molecular_mass:
                molecular_mass = response_dict.get("molecularMass", None)
                if molecular_mass:
                    try:
                        self.molecular_mass = float(molecular_mass)
                    except ValueError as e:
                        logger.warn(
                            f"Could not convert molecular mass results'{molecular_mass}'"
                            + "returned from CAS api request.",
                            exc_info=e)

            experimental_properties = response_dict.get("experimentalProperties", [])
            props = [
                CASExperimentalProperty.m_from_dict(p) for p in experimental_properties]
            if len(props) > 0:
                self.cas_experimental_properties = props

            property_citations = response_dict.get("propertyCitations", [])
            citations = [CASPropertyCitation.m_from_dict(c) for c in property_citations]
            if len(citations) > 0:
                self.cas_property_citations = citations

            self.cas_synonyms = response_dict.get("synonyms", [])

        elif response.status_code == 404:
            logger.warn(f"No CAS entry found with CAS number: {self.cas_number}")
        elif response.status_code >= 500:
            logger.warn(f"Remote server error on CAS API call.")
        else:
            logger.warn(
                f"Unexpected response code: {response.status_code} from CAS API call.")

    def _cas_search_unique(self, search: str, archive, logger: Any) -> bool:
        '''Private method for performing a search of the CAS API and populating the
        attributes with the CAS number of any unique search result.

        Args:
            search (str): The string to search the CAS API with.
            archive (EntryArchive): The archive that is being normalized.
            logger (Any): A structlog logger.

        Returns:
            bool: Whether the search found a unique result.
        '''
        import httpx
        response = httpx.get(
            f'https://commonchemistry.cas.org/api/search?q={search}', timeout=2)
        if response.status_code == 200:
            response_dict = response.json()
            n_hits = response_dict.get("count", 0)
            if n_hits == 0:
                logger.info(f"Queried the CAS API for '{search}' without result.")
                return False
            elif n_hits == 1:
                try:
                    self.cas_number = response_dict["results"][0]["rn"]
                    self._populate_from_cas(archive, logger)
                    return True
                except KeyError as e:
                    logger.warn("Unknown response format from CAS API.", exc_info=e)
            else:
                logger.warn(
                    f"Query to CAS API for '{search}' returned {n_hits} hits, please "
                    + "specify further.")
                return False
        return False

    def _find_cas(self, archive, logger: Any) -> None:
        '''Private method for finding the CAS number using the filled attributes in the
        following order:
        1. `cas_name`
        2. `inchi`
        3. `inchi_key`
        4. `smile`
        5. `canonical_smile`
        6. `name`
        The first unique hit will populate the `cas_number` attribute and return.

        Args:
            archive (EntryArchive): The archive that is being normalized.
            logger (Any): A structlog logger.
        '''
        for key in (
                self.cas_name,
                self.inchi,
                self.inchi_key,
                self.smile,
                self.canonical_smile,
                self.name):
            if key and self._cas_search_unique(key, archive, logger):
                return

    def normalize(self, archive, logger: Any) -> None:
        '''The normalizer method for the `Substance` class.
        This method will attempt to get data on the substance instance from the CAS API:
        https://commonchemistry.cas.org/api-overview
        If a CAS number is specified the details are retrieved directly.
        Otherwise a search query is made for the filled attributes in the following order:
        1. `cas_name`
        2. `inchi`
        3. `inchi_key`
        4. `smile`
        5. `canonical_smile`
        6. `name`

        Args:
            archive (EntryArchive): The archive that is being normalized.
            logger (Any): A structlog logger.
        '''
        super(Substance, self).normalize(archive, logger)
        if logger is None:
            logger = utils.get_logger(__name__)
        from nomad.atomutils import Formula

        if self.cas_number:
            self.cas_number = self.cas_number.strip()
            self._populate_from_cas(archive, logger)
        else:
            self._find_cas(archive, logger)

        if self.molecular_formula:
            if not archive.results:
                archive.results = Results()
            if not archive.results.material:
                archive.results.material = Material()

            try:
                formula = Formula(self.molecular_formula)
                formula.populate_material(material=archive.results.material)
                if not self.elemental_composition:
                    for element, fraction in formula.atomic_fractions().items():
                        self.elemental_composition.append(ElementalComposition.m_from_dict({
                            "element": element, "atomic_fraction": fraction}))
            except Exception as e:
                logger.warn('Could not analyse chemical formula.', exc_info=e)

        super(Substance, self).normalize(archive, logger)


class SampleID(ArchiveSection):
    ''' A ELN base section that can be used for sample IDs.
    If the `sample_owner`, `sample_short_name`, `ìnstitute`, and `creation_datetime`
    quantities are provided, the sample_id will be automatically created as a combination
    of these four quantities.'''

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

    def normalize(self, archive, logger):
        '''The normalizer for the `SampleID` class.
        If sample owner is not filled the field will be filled by the first two letters of
        the first name joined with the first two letters of the last name of the author.
        If the institute is not filled a institute abreviations will be constructed from
        the author's affiliation.
        If no creation datetime is filled, the datetime will be taken from the `datetime`
        property of the parent, if it exists, otherwise the current date and time will be
        used.
        If no short name is filled, the name will be taken from the parent name, if it
        exists, otherwise it will be taken from the archive metadata entry name, if it
        exists, and finally if no other options are available it will use the name of the
        mainfile.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (Any): A structlog logger.
        '''
        super(SampleID, self).normalize(archive, logger)

        if self.sample_owner is None or self.institute is None:
            from unidecode import unidecode
            author = archive.metadata.main_author
            if author and self.sample_owner is None:
                first_short = unidecode(author.first_name)[:2]
                last_short = unidecode(author.last_name)[:2]
                self.sample_owner = first_short + last_short
            if author and self.institute is None:
                unwanted_words = ("zu", "of", "the", "fur", "für")
                institute = ''
                all_words = re.split(' |-|_|,|:|;', unidecode(author.affiliation))
                wanted_words = [w for w in all_words if w.lower() not in unwanted_words]
                for word in wanted_words:
                    word_remainder = word
                    while word_remainder and len(institute) < 3:
                        letter = word_remainder[:1]
                        if letter.isupper():
                            institute += letter
                        word_remainder = word_remainder[1:]
                if len(institute) < min(len(author.affiliation), 2):
                    if len(wanted_words) < 3:
                        institute = author.affiliation[:3].upper()
                    else:
                        institute = ''.join([w[:1] for w in wanted_words[:3]]).upper()
                self.institute = institute

        if self.creation_datetime is None:
            if self.m_parent and getattr(self.m_parent, 'datetime', None):
                self.creation_datetime = self.m_parent.datetime
            else:
                self.creation_datetime = datetime.datetime.now()

        if self.sample_short_name is None:
            if self.m_parent and getattr(self.m_parent, 'name', None):
                name = self.m_parent.name
            elif archive.metadata.entry_name:
                name = archive.metadata.entry_name
            else:
                name = archive.metadata.mainfile
            self.sample_short_name = re.sub(r'_|\s', '-', name.split('.')[0])

        if self.institute and self.sample_short_name and self.sample_owner and self.creation_datetime:
            creation_date = self.creation_datetime.strftime('%Y%m%d')
            sample_owner = self.sample_owner.replace(' ', '-')
            sample_id_list = [self.institute, sample_owner, creation_date, self.sample_short_name]
            self.sample_id = '_'.join(sample_id_list)

        if not archive.results:
            archive.results = Results(eln=ELN())
        if not archive.results.eln:
            archive.results.eln = ELN()

        if self.sample_id:
            if not archive.results.eln.lab_ids:
                archive.results.eln.lab_ids = []
            if self.sample_id not in archive.results.eln.lab_ids:
                archive.results.eln.lab_ids.append(self.sample_id)

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
        import httpx

        # Parse journal name, lead author and publication date from crossref
        if self.DOI_number:
            try:
                url = f'https://api.crossref.org/works/{self.DOI_number}?mailto=contact@nomad-lab.eu'
                timeout = 2
                r = httpx.get(url, timeout=timeout)
                if r.status_code == 200:
                    temp_dict = r.json()
                    # make sure the doi has the prefix https://doi.org/
                    if self.DOI_number.startswith('10.'):
                        self.DOI_number = 'https://doi.org/' + self.DOI_number
                    self.publication_authors = [f"{v['given']} {v['family']}" for v in temp_dict['message']['author']]
                    self.journal = temp_dict['message']['container-title'][0]
                    self.publication_title = temp_dict['message']['title'][0]
                    self.publication_date = dateutil.parser.parse(temp_dict['message']['created']['date-time'])
                    if not archive.metadata:
                        archive.metadata = EntryMetadata()
                    if not archive.metadata.references:
                        archive.metadata.references = []
                    # if any item in the references list starts with 10. add the prefix https://doi.org/
                    for i, ref in enumerate(archive.metadata.references):
                        if ref.startswith('10.'):
                            archive.metadata.references[i] = 'https://doi.org/' + ref
                    if self.DOI_number not in archive.metadata.references:
                        archive.metadata.references.append(self.DOI_number)

                else:
                    logger.warning(f'Could not parse DOI number {self.DOI_number}')
            except Exception as e:
                logger.warning(f'Could not parse crossref for {self.DOI_number}')
                logger.warning(e)


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
                formula.populate_material(material=archive.results.material)
            except Exception as e:
                logger.warn('could not analyse chemical formula', exc_info=e)


class Chemical(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for chemicals.'''
    pass


class Sample(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for samples.'''
    pass


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
