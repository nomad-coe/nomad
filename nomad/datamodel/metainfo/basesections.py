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

import datetime
import re
from typing import (
    Dict,
    List,
)

import numpy as np
from ase.data import (
    chemical_symbols,
    atomic_numbers,
    atomic_masses,
)
from structlog.stdlib import (
    BoundLogger,
)

from nomad import (
    utils,
)
from nomad.units import (
    ureg,
)
from nomad.metainfo import (
    Quantity,
    Datetime,
    Reference,
    Section,
    SubSection,
)
from nomad.metainfo.util import (
    MEnum,
)
from nomad.datamodel.data import (
    ArchiveSection,
    EntryData,
)
from nomad.datamodel.results import (
    Results,
    ELN,
    ElementalComposition as ResultsElementalComposition,
    Material,
)


class BaseSection(ArchiveSection):
    '''
    A generic abstract base section that provides a few commonly used properties.

    If you inherit from this section, but do not need some quantities, list those
    quantities in the `eln.hide` annotation of your inheriting section definition.

    Besides predefining some quantities, these base sections will add some metadata
    to NOMAD's search. A particular example are `tags`, if you define a string
    or enum quantity in your sections named `tags`, its values will be searchable.
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/BFO_0000001'],
    )
    name = Quantity(
        type=str,
        description='A short human readable and descriptive name.',
        a_eln=dict(component='StringEditQuantity', label='Short name'),
    )
    datetime = Quantity(
        type=Datetime,
        description='The date and time associated with this section.',
        a_eln=dict(component='DateTimeEditQuantity'),
    )
    lab_id = Quantity(
        type=str,
        description='''An ID string that is unique at least for the lab that produced this
            data.''',
        a_eln=dict(component='StringEditQuantity', label="ID"),
    )
    description = Quantity(
        type=str,
        description='Any information that cannot be captured in the other fields.',
        a_eln=dict(component='RichTextEditQuantity'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `BaseSection` class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(BaseSection, self).normalize(archive, logger)

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


class Entity(BaseSection):
    '''
    An object that persists, endures, or continues to exist through time while maintaining
    its identity.
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/BFO_0000002'],
    )


class Activity(BaseSection):
    '''
    An action that has a temporal extension and for some time depends on some entity.
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/BFO_0000015'],
    )
    datetime = Quantity(
        type=Datetime,
        description='The date and time when this activity was started.',
        a_eln=dict(component='DateTimeEditQuantity', label='Starting Time'),
    )
    end_time = Quantity(
        type=Datetime,
        description='The date and time when this activity was done.',
        a_eln=dict(component='DateTimeEditQuantity', label='Ending Time'),
    )
    method = Quantity(
        type=str,
        description='A short consistent handle for the applied method.',
    )
    location = Quantity(
        type=str,
        description='location of the experiment.',
        a_eln=dict(component='StringEditQuantity'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `Activity` class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(Activity, self).normalize(archive, logger)

        if archive.results.eln.methods is None:
            archive.results.eln.methods = []
        if self.method:
            archive.results.eln.methods.append(self.method)
        else:
            archive.results.eln.methods.append(self.m_def.name)


class ElementalComposition(ArchiveSection):
    '''
    A section for describing the elemental composition of a system, i.e. the element
    and its atomic fraction.
    '''
    m_def = Section(
        label_quantity='element',
    )
    element = Quantity(
        type=MEnum(chemical_symbols[1:]),
        description='''
        The symbol of the element, e.g. 'Pb'.
        ''',
        a_eln=dict(component='AutocompleteEditQuantity'),
    )
    atomic_fraction = Quantity(
        type=np.float64,
        description='''
        The atomic fraction of the element in the system it is contained within.
        Per definition a positive value less than or equal to 1.
        ''',
        a_eln=dict(component='NumberEditQuantity'),
    )
    mass_fraction = Quantity(
        type=np.float64,
        description='''
        The mass fraction of the element in the system it is contained within.
        Per definition a positive value less than or equal to 1.
        ''',
        a_eln=dict(component='NumberEditQuantity'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `ElementalComposition` class. Will add a
        results.material subsection if none exists. Will append the element to the
        elements property of that subsection and a
        nomad.datamodel.results.ElementalComposition instances to the
        elemental_composition property  using the element and atomic fraction from this
        section.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(ElementalComposition, self).normalize(archive, logger)

        if self.element:
            if not archive.results:
                archive.results = Results()
            if not archive.results.material:
                archive.results.material = Material()

            if self.element not in chemical_symbols:
                logger.warn(
                    f"'{self.element}' is not a valid element symbol and this "
                    "elemental_composition section will be ignored.")
            elif self.element not in archive.results.material.elements:
                archive.results.material.elements += [self.element]
            if self.atomic_fraction or self.mass_fraction:
                comp_result_section = archive.results.material.elemental_composition
                result_composition = ResultsElementalComposition(
                    element=self.element,
                    atomic_fraction=self.atomic_fraction,
                    mass_fraction=self.mass_fraction,
                    mass=atomic_masses[atomic_numbers[self.element]] * ureg.amu,
                )
                existing_elements = [comp.element for comp in comp_result_section]
                if self.element in existing_elements:
                    index = existing_elements.index(self.element)
                    comp_result_section[index].atomic_fraction = self.atomic_fraction
                    comp_result_section[index].mass_fraction = self.mass_fraction
                    comp_result_section[index].mass = atomic_masses[atomic_numbers[self.element]] * ureg.amu
                else:
                    comp_result_section.append(result_composition)


class System(Entity):
    '''
    A base section for any system of materials which is investigated or used to construct
    other systems.
    '''
    elemental_composition = SubSection(
        description='''
        A list of all the elements found in the system together and their respective
        atomic fraction within the system.
        ''',
        section_def=ElementalComposition,
        repeats=True,
    )

    def _fill_fractions(self, archive, logger: BoundLogger) -> None:
        '''
        Private method for attempting to fill missing fractions (atomic or mass) in the
        `ElementalComposition` objects listed in elemental_composition.
        If a single fraction (of mass or atomic) is left blank it will be calculated from
        the others. If after this check, all atomic fractions are filled and no mass
        fractions are, the mass fractions will be calculated from the atomic fractions.
        similarly, the atomic fractions are calculated from the mass ones if all mass and
        no atomic fractions are filled.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.

        Raises:
            ValueError: If an unknown element is present.
        '''
        atomic_fractions = []
        mass_fractions = []
        element_masses = []
        for comp in self.elemental_composition:
            atomic_fractions.append(comp.atomic_fraction)
            mass_fractions.append(comp.mass_fraction)
            if comp.element is None:
                return
            elif comp.element not in atomic_numbers:
                raise ValueError(f'Unknown element symbol: {comp.element}')
            element_masses.append(atomic_masses[atomic_numbers[comp.element]])
        atomic_blanks = atomic_fractions.count(None)
        mass_blanks = mass_fractions.count(None)
        if atomic_blanks == 1:
            index = atomic_fractions.index(None)
            balance = 1 - np.nansum(np.array(atomic_fractions, dtype=np.float64))
            self.elemental_composition[index].atomic_fraction = balance
            self.elemental_composition[index].normalize(archive, logger)
            atomic_fractions[index] = balance
            atomic_blanks = 0
        if mass_blanks == 1:
            index = mass_fractions.index(None)
            balance = 1 - np.nansum(np.array(mass_fractions, dtype=np.float64))
            self.elemental_composition[index].mass_fraction = balance
            self.elemental_composition[index].normalize(archive, logger)
            mass_fractions[index] = balance
            mass_blanks = 0
        if atomic_blanks == 0 and mass_blanks == len(mass_fractions):
            masses = np.array(atomic_fractions) * np.array(element_masses)
            mass_fractions = masses / masses.sum()
            for index, mass_fraction in enumerate(mass_fractions):
                self.elemental_composition[index].mass_fraction = mass_fraction
                self.elemental_composition[index].normalize(archive, logger)
        if mass_blanks == 0 and atomic_blanks == len(atomic_fractions):
            atoms = np.array(mass_fractions) / np.array(element_masses)
            atomic_fractions = atoms / atoms.sum()
            for index, atomic_fraction in enumerate(atomic_fractions):
                self.elemental_composition[index].atomic_fraction = atomic_fraction
                self.elemental_composition[index].normalize(archive, logger)

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `System` class. Will attempt to fill mass fractions or
        atomic fractions if left blank.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(System, self).normalize(archive, logger)

        if len(self.elemental_composition) > 0:
            self._fill_fractions(archive, logger)


class Instrument(Entity):
    '''
    A base section that can be used for instruments.
    '''
    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `Instrument` class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(Instrument, self).normalize(archive, logger)

        if self.name:
            if archive.results.eln.instruments is None:
                archive.results.eln.instruments = []
            archive.results.eln.instruments.append(self.name)


class Process(Activity):
    '''
    A planned process which results in physical changes in a specified input material.
    [ obi : prs obi : mc obi : fg obi : jf obi : bp ]

    Synonyms:
     - preparative method
     - sample preparation
     - sample preparative method
     - material transformations
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/OBI_0000094'],
    )


class Measurement(Activity):
    '''
    A planned process with the objective to produce information about the material entity
    that is the evaluant, by physically examining it or its proxies. [ obi : pppb ]
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/OBI_0000070'],
    )


class SynthesisMethod(Process):
    '''
    A method used to synthesise a sample.
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/CHMO_0001301'],
    )


class Component(ArchiveSection):
    '''
    A section for describing a component and its role in an ensemble.
    '''
    name = Quantity(
        type=str,
        description='A short name for the component.',
        a_eln=dict(component='StringEditQuantity', label='Component name'),
    )
    system = Quantity(
        type=Reference(System.m_def),
        description='A reference to the component system.',
        a_eln=dict(component='ReferenceEditQuantity'),
    )
    mass = Quantity(
        type=np.float64,
        description='The mass of the component.',
        unit='kg',
        a_eln=dict(component='NumberEditQuantity', defaultDisplayUnit='mg'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `Component` class. If none is set, the normalizer
        will set the name of the component to be that of the referenced system if it has
        one.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(Component, self).normalize(archive, logger)
        if self.name is None and self.system is not None:
            self.name = self.system.name


class Ensemble(System):
    '''
    A base section for an ensemble of material systems. Each component of the ensemble is
    of a (sub)type of `System`.
    '''
    components = SubSection(
        description='''
        A list of all the components of the ensemble containing a name, reference to the
        system section and mass of that component.
        ''',
        section_def=Component,
        repeats=True,
    )

    @staticmethod
    def _atomic_to_mass(composition: List[ElementalComposition], mass: float) -> Dict[str, float]:
        '''
        Private static method for converting list of ElementalComposition objects to
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
        '''
        Private static method for converting dictionary of elements with their masses
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

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `Ensemble` class. If the elemental composition list is
        empty, the normalizer will iterate over the components and extract all the
        elements for populating the elemental composition list. If masses are provided for
        all components and the elemental composition of all components contain atomic
        fractions the normalizer will also calculate the atomic fractions for the
        ensemble. The populated elemental composition list is added to the results by the
        normalizer in the `System` super class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
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
        for comp in self.elemental_composition:
            comp.normalize(archive, logger)

        super(Ensemble, self).normalize(archive, logger)


class Substance(System):
    '''
    A base section for any elemental, molecular of single phase substance.
    '''
    m_def = Section(
        links=['http://purl.obolibrary.org/obo/CHEBI_23367'],
    )
    name = Quantity(
        type=str,
        description='The name of the substance entry.',
        a_eln=dict(component='StringEditQuantity', label='Substance name'),
    )
    lab_id = Quantity(
        type=str,
        description='''
        A human human readable substance ID that is at least unique for the lab.
        ''',
        a_eln=dict(component='StringEditQuantity', label='Substance ID'),
    )
    molecular_formula = Quantity(
        type=str,
        description='Molecular formula.',
        a_eln=dict(component='StringEditQuantity'),
    )
    description = Quantity(
        type=str,
        description='''
        A field for adding additional information about the substance that is not captured
        by the other quantities and subsections.
        ''',
        a_eln=dict(
            component='RichTextEditQuantity',
            label='Detailed substance description',
        ),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer method for the `Substance` class.
        This method will populate the results.material section and the elemental
        composition sub section using the molecular formula.

        Args:
            archive (EntryArchive): The archive that is being normalized.
            logger (BoundLogger): A structlog logger.
        '''
        super(Substance, self).normalize(archive, logger)
        if logger is None:
            logger = utils.get_logger(__name__)
        from nomad.atomutils import Formula

        if self.molecular_formula:
            if not archive.results:
                archive.results = Results()
            if not archive.results.material:
                archive.results.material = Material()

            try:
                formula = Formula(self.molecular_formula)
                formula.populate(archive.results.material)
                if not self.elemental_composition:
                    mass_fractions = formula.mass_fractions()
                    for element, fraction in formula.atomic_fractions().items():
                        self.elemental_composition.append(
                            ElementalComposition(
                                element=element,
                                atomic_fraction=fraction,
                                mass_fraction=mass_fractions[element],
                            )
                        )
            except Exception as e:
                logger.warn('Could not analyse chemical formula.', exc_info=e)


class CASExperimentalProperty(ArchiveSection):
    '''
    A section for experimental properties retrieved from the CAS API.
    '''
    name = Quantity(
        type=str,
        description='CAS experimental property name.',
    )
    property = Quantity(
        type=str,
        description='CAS experimental property.',
    )
    sourceNumber = Quantity(
        type=str,
        description='CAS experimental property source.',
    )


class CASPropertyCitation(ArchiveSection):
    '''
    A section for citations of the experimental properties retrieved from the CAS API.
    '''
    docUri = Quantity(
        type=str,
        description='CAS property citation document uri.',
    )
    sourceNumber = Quantity(
        type=int,
        decription='CAS property citation source number.',
    )
    source = Quantity(
        type=str,
        description='CAS property citation source.',
    )


class CASSubstance(Substance):
    '''
    A base section for any Substance with a CAS number.
    '''
    cas_uri = Quantity(
        type=str,
        description='CAS uri',
        a_eln=dict(component='StringEditQuantity', label='CAS uri'),
    )
    cas_number = Quantity(
        type=str,
        description='CAS number.',
        a_eln=dict(component='StringEditQuantity', label='CAS number'),
    )
    cas_name = Quantity(
        type=str,
        description='CAS name.',
        a_eln=dict(component='StringEditQuantity', label='CAS name'),
    )
    image = Quantity(
        type=str,
        description='CAS image.',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor', label='Image of substance'),
    )
    inchi = Quantity(
        type=str,
        description='CAS inchi.',
        a_eln=dict(component='StringEditQuantity'),
    )
    inchi_key = Quantity(
        type=str,
        description='CAS inchi key.',
        a_eln=dict(component='StringEditQuantity'),
    )
    smile = Quantity(
        type=str,
        description='CAS smile.',
        a_eln=dict(component='StringEditQuantity'),
    )
    canonical_smile = Quantity(
        type=str,
        description='CAS canonical smile.',
        a_eln=dict(component='StringEditQuantity'),
    )
    molecular_formula = Quantity(
        type=str,
        description='CAS molecular formula.',
        a_eln=dict(component='StringEditQuantity'),
    )
    molecular_mass = Quantity(
        type=np.dtype(np.float64),
        unit='Da',
        description='CAS molecular mass.',
        a_eln=dict(component='NumberEditQuantity'),
    )
    cas_experimental_properties = SubSection(
        section_def=CASExperimentalProperty,
        repeats=True,
    )
    cas_property_citations = SubSection(
        section_def=CASPropertyCitation,
        repeats=True,
    )
    cas_synonyms = Quantity(
        type=str,
        shape=['*'],
        description='CAS synonyms.',
    )

    def _populate_from_cas(self, archive, logger: BoundLogger) -> None:
        '''
        Private method for populating the attributes from a call to the CAS API using
        the `cas_number`.
        Will overwrite exisiting CAS attributes if the query provides a value for them.
        I.e. all attributes that begin with `cas_`.

        Args:
            archive (EntryArchive): The archive that is being normalized.
            logger (BoundLogger): A structlog logger.
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

    def _cas_search_unique(self, search: str, archive, logger: BoundLogger) -> bool:
        '''
        Private method for performing a search of the CAS API and populating the
        attributes with the CAS number of any unique search result.

        Args:
            search (str): The string to search the CAS API with.
            archive (EntryArchive): The archive that is being normalized.
            logger (BoundLogger): A structlog logger.

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

    def _find_cas(self, archive, logger: BoundLogger) -> None:
        '''
        Private method for finding the CAS number using the filled attributes in the
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
            logger (BoundLogger): A structlog logger.
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

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer method for the `CASSubstance` class.
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
            logger (BoundLogger): A structlog logger.
        '''
        super(CASSubstance, self).normalize(archive, logger)
        if logger is None:
            logger = utils.get_logger(__name__)

        if self.cas_number:
            self.cas_number = self.cas_number.strip()
            self._populate_from_cas(archive, logger)
        else:
            self._find_cas(archive, logger)

        super(CASSubstance, self).normalize(archive, logger)


class SampleID(ArchiveSection):
    '''
    A base section that can be used for sample IDs.
    If the `sample_owner`, `sample_short_name`, `ìnstitute`, and `creation_datetime`
    quantities are provided, the sample_id will be automatically created as a combination
    of these four quantities.
    '''
    institute = Quantity(
        type=str,
        description='Alias/short name of the home institute of the owner, i.e. *HZB*.',
        a_eln=dict(component='StringEditQuantity'),
    )
    sample_owner = Quantity(
        type=str,
        shape=[],
        description='Name or alias of the process operator, e.g. jmp',
        a_eln=dict(component='StringEditQuantity'),
    )
    creation_datetime = Quantity(
        type=Datetime,
        description='Creation date of the sample.',
        a_eln=dict(component='DateTimeEditQuantity'),
    )
    sample_short_name = Quantity(
        type=str,
        description='''A short name of the sample (the identifier scribed on the smaple,
         or in the sample container), e.g. 4001-8, YAG-2-34.
         This is to be managed and decided internally by the labs,
         although we recomend to avoid the following characters on it: "_", "/", "\" and "."''',
        a_eln=dict(component='StringEditQuantity'),
    )
    sample_id = Quantity(
        type=str,
        description='''Full sample id. Ideally a human readable sample id convention,
        which is simple, understandable and still having chances of becoming unique.
        If the `sample_owner`, `sample_short_name`, `ìnstitute`, and `creation_datetime`
        are provided, this will be formed automatically by joining these components by an underscore (_).
        Spaces in any of the individual components will be replaced with hyphens (-).
        An example would be hzb_oah_20200602_4001-08''',
        a_eln=dict(component='StringEditQuantity'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `SampleID` class.
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
            logger (BoundLogger): A structlog logger.
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
    '''
    A base section that can be used for references.
    '''
    DOI_number = Quantity(
        type=str,
        shape=[],
        description="""
            The DOI number referring to the published paper or dataset where the data can be found.
            Examples:
            10.1021/jp5126624
            10.1016/j.electacta.2017.06.032
        """,
        a_eln=dict(component='EnumEditQuantity', props=dict(suggestions=[])),
    )

    publication_authors = Quantity(
        type=str,
        shape=['*'],
        description="""
            The authors of the publication.
            If several authors, end with et al. If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """,
    )
    publication_date = Quantity(
        type=Datetime,
        shape=[],
        description="""
            Publication date.
            If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """,
    )
    journal = Quantity(
        type=str,
        shape=[],
        description="""
            Name of the journal where the data is published.
            If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """,
    )
    publication_title = Quantity(
        type=str,
        shape=[],
        description="""
            Title of the publication.
            If the DOI number is given correctly,
            this will be extracted automatically from www.crossref.org
        """,
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        '''
        The normalizer for the `PublicationReference` class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        '''
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
                logger.warning(str(e))
