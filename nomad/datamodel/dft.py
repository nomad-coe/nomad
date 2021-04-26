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

'''
DFT specific metadata
'''

import re

from nomad import config, utils
from nomad.metainfo import MSection, Section, Quantity, MEnum, SubSection
from nomad.metainfo.search_extension import Search

from .optimade import OptimadeEntry
from .metainfo.common_dft import Workflow, FastAccess
from .metainfo.common_dft import section_XC_functionals


xc_treatments = {
    'gga': 'GGA',
    'hf_': 'HF',
    'oep': 'OEP',
    'hyb': 'hybrid',
    'mgg': 'meta-GGA',
    'vdw': 'vdW',
    'lda': 'LDA',
}
''' https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional '''

basis_sets = {
    'gaussians': 'gaussians',
    'realspacegrid': 'real-space grid',
    'planewaves': 'plane waves'
}

compound_types = [
    'unary',
    'binary',
    'ternary',
    'quaternary',
    'quinary',
    'sexinary',
    'septenary',
    'octanary',
    'nonary',
    'decinary'
]

_electronic_quantities = [
    'electronic_band_structure',
    'electronic_dos',
    'eigenvalues_values',
]

_mechanical_quantities = [
    'stress_tensor'
]

_thermal_quantities = [
    'thermodynamical_property_heat_capacity_C_v',
    'vibrational_free_energy_at_constant_volume',
    'phonon_band_structure',
    'phonon_dos',
]

_magnetic_quantities = [
    'spin_S2'
]

_optical_quantities = [
    'oscillator_strengths',
    'transition_dipole_moments'
]

_searchable_quantities = set(_electronic_quantities + _mechanical_quantities + _thermal_quantities + _magnetic_quantities + _optical_quantities)

version_re = re.compile(r'(\d+(\.\d+(\.\d+)?)?)')


def map_functional_name_to_xc_treatment(name):
    if name == config.services.unavailable_value:
        return name

    return xc_treatments.get(name[:3].lower(), config.services.unavailable_value)


def map_basis_set_to_basis_set_label(name):
    key = name.replace('_', '').replace('-', '').replace(' ', '').lower()
    return basis_sets.get(key, name)


def simplify_version(version):
    match = version_re.search(version)
    if match is None:
        return version
    else:
        return match.group(0)


def valid_array(array):
    """Checks if the given variable is a non-empty array.
    """
    return array is not None and len(array) > 0


class Label(MSection):
    '''
    Label that further classify a structure.

    Attributes:
        label: The label as a string
        type: The type of the label
        source: The source that this label was taken from.

    '''
    label = Quantity(type=str, a_search=Search())

    type = Quantity(type=MEnum(
        'compound_class', 'classification', 'prototype', 'prototype_id'),
        a_search=Search())

    source = Quantity(
        type=MEnum('springer', 'aflow_prototype_library'),
        a_search=Search())


class DFTMetadata(MSection):
    m_def = Section(a_domain='dft')

    basis_set = Quantity(
        type=str, default='not processed',
        description='The used basis set functions.',
        a_search=Search(statistic_values=[
            '(L)APW+lo', 'gaussians', 'numeric AOs', 'plane waves', 'psinc functions',
            'real-space grid', 'unavailable', 'not processed'
        ]))

    xc_functional = Quantity(
        type=str, default='not processed',
        description='The libXC based xc functional classification used in the simulation.',
        a_search=Search(
            statistic_values=list(xc_treatments.values()) + ['unavailable', 'not processed'],
            statistic_size=100))

    xc_functional_names = Quantity(
        type=str, default=[], shape=['*'],
        description='The list of libXC functional names that where used in this entry.',
        a_search=Search(many_and='append'))

    system = Quantity(
        type=str, default='not processed',
        description='The system type of the simulated system.',
        a_search=Search(statistic_values=[
            '1D', '2D', 'atom', 'bulk', 'molecule / cluster', 'surface',
            'unavailable', 'not processed'
        ]))

    compound_type = Quantity(
        type=str, default='not processed',
        description='The compound type of the simulated system.',
        a_search=Search(statistic_values=compound_types + ['not processed'])
    )

    crystal_system = Quantity(
        type=str, default='not processed',
        description='The crystal system type of the simulated system.',
        a_search=Search(
            statistic_values=[
                'cubic', 'hexagonal', 'monoclinic', 'orthorhombic', 'tetragonal',
                'triclinic', 'trigonal', 'unavailable', 'not processed']
        ))

    spacegroup = Quantity(
        type=int, default=-1,
        description='The spacegroup of the simulated system as number.',
        a_search=Search())

    spacegroup_symbol = Quantity(
        type=str, default='not processed',
        description='The spacegroup as international short symbol.',
        a_search=Search())

    code_name = Quantity(
        type=str, default='not processed',
        description='The name of the used code.',
        a_search=Search())  # in import the parser module is added codes here as statistic_values

    code_version = Quantity(
        type=str, default='not processed',
        description='The version of the used code.',
        a_search=Search())

    n_geometries = Quantity(
        type=int, default=0, description='Number of unique geometries.',
        a_search=Search(metric_name='geometries', metric='sum'))

    n_calculations = Quantity(
        type=int, default=0,
        description='Number of single configuration calculation sections',
        a_search=Search(metric_name='calculations', metric='sum'))

    n_total_energies = Quantity(
        type=int, default=0, description='Number of total energy calculations',
        a_search=Search(metric_name='total_energies', metric='sum'))

    n_quantities = Quantity(
        type=int, default=0, description='Number of metainfo quantities parsed from the entry.',
        a_search=Search(metric='sum', metric_name='quantities'))

    quantities = Quantity(
        type=str, shape=['0..*'],
        description='All quantities that are used by this entry.',
        a_search=Search(
            metric_name='distinct_quantities', metric='cardinality', many_and='append'))

    searchable_quantities = Quantity(
        type=str, shape=['0..*'],
        description='All quantities with existence filters in the search GUI.',
        a_search=Search(many_and='append', statistic_size=len(_searchable_quantities) + 15))  # Temporarily increased the statistics size while migrating from old set to new one.

    geometries = Quantity(
        type=str, shape=['0..*'],
        description='Hashes for each simulated geometry',
        a_search=Search(metric_name='unique_geometries', metric='cardinality'))

    group_hash = Quantity(
        type=str,
        description='Hashes that describe unique geometries simulated by this code run.',
        a_search=Search(many_or='append', group='groups_grouped', metric_name='groups', metric='cardinality'))

    labels = SubSection(
        sub_section=Label, repeats=True, categories=[FastAccess],
        description='The labels taken from AFLOW prototypes and springer.',
        a_search=Search())

    labels_springer_compound_class = Quantity(
        type=str, shape=['0..*'],
        description='Springer compund classification.',
        a_search=Search(
            many_and='append', statistic_size=20,
            statistic_order='_count'))

    labels_springer_classification = Quantity(
        type=str, shape=['0..*'],
        description='Springer classification by property.',
        a_search=Search(
            many_and='append', statistic_size=20,
            statistic_order='_count'))

    optimade = SubSection(
        sub_section=OptimadeEntry,
        description='Metadata used for the optimade API.',
        a_search=Search())

    workflow = Quantity(type=Workflow, a_search=Search())

    def code_name_from_parser(self):
        entry = self.m_parent
        if entry.parser_name is not None:
            from nomad.parsing.parsers import parser_dict
            parser = parser_dict.get(entry.parser_name)
            if hasattr(parser, 'code_name'):
                return parser.code_name
        return config.services.unavailable_value

    def update_group_hash(self):
        user_id = None
        uploader = self.m_parent.uploader
        if uploader is not None:
            user_id = uploader.user_id
        self.group_hash = utils.hash(
            self.m_parent.formula,
            self.spacegroup,
            self.basis_set,
            self.xc_functional,
            self.code_name,
            self.code_version,
            self.m_parent.with_embargo,
            user_id)

    def apply_domain_metadata(self, entry_archive):
        from nomad.normalizing.system import normalized_atom_labels
        entry = self.m_parent

        logger = utils.get_logger(__name__).bind(
            upload_id=entry.upload_id, calc_id=entry.calc_id, mainfile=entry.mainfile)

        self.code_name = self.code_name_from_parser()

        if entry_archive is None:
            return

        section_run = entry_archive.section_run
        if not section_run:
            logger.warn('no section_run found')
            return
        section_run = section_run[0]

        # default values
        self.system = config.services.unavailable_value
        self.crystal_system = config.services.unavailable_value
        self.spacegroup_symbol = config.services.unavailable_value
        self.basis_set = config.services.unavailable_value
        self.xc_functional = config.services.unavailable_value

        section_system = None
        for section in section_run.section_system:
            if section.is_representative:
                section_system = section
                break

        # code and code specific ids
        try:
            code_name = section_run.program_name
            if code_name:
                self.code_name = code_name
            else:
                raise KeyError
        except KeyError as e:
            logger.warn('archive without program_name', exc_info=e)

        try:
            version = section_run.program_version
            if version:
                self.code_version = simplify_version(version)
            else:
                raise KeyError
        except KeyError:
            self.code_version = config.services.unavailable_value

        def get_value(value):
            return value if value else config.services.unavailable_value

        raw_id = section_run.raw_id
        if raw_id is not None:
            entry.raw_id = raw_id

        # metadata (system, method, chemistry)
        atom_labels = section_system.atom_labels if section_system else []
        atoms = atom_labels if atom_labels else []
        entry.n_atoms = len(atoms)
        atoms = list(set(normalized_atom_labels(set(atoms))))
        atoms.sort()
        entry.atoms = atoms
        self.compound_type = compound_types[len(atoms) - 1] if len(atoms) <= 10 else '>decinary'

        self.system = config.services.unavailable_value
        self.crystal_system = config.services.unavailable_value
        self.spacegroup_symbol = config.services.unavailable_value

        section_symmetry = None
        if section_system and len(section_system.section_symmetry) > 0:
            section_symmetry = section_system.section_symmetry[0]
            self.crystal_system = get_value(section_symmetry.crystal_system)
            spacegroup = section_symmetry.space_group_number
            self.spacegroup = 0 if not spacegroup else int(spacegroup)
            self.spacegroup_symbol = get_value(section_symmetry.international_short_symbol)

        program_basis_set_type = section_run.program_basis_set_type
        if program_basis_set_type:
            self.basis_set = map_basis_set_to_basis_set_label(program_basis_set_type)

        if section_system:
            self.system = get_value(section_system.system_type)
            entry.formula = get_value(section_system.chemical_composition_bulk_reduced)

        # metrics and quantities
        quantities = set()
        searchable_quantities = set()
        geometries = set()
        xc_functionals = set()
        xc_functional = None

        n_quantities = 0
        n_calculations = 0
        n_total_energies = 0
        n_geometries = 0

        for section, property_def, _ in entry_archive.m_traverse():
            property_name = property_def.name
            quantities.add(property_name)
            n_quantities += 1

            if property_name in _searchable_quantities:
                searchable_quantities.add(property_name)

            if property_def == section_XC_functionals.XC_functional_name:
                xc_functional = getattr(section, property_name)
                if xc_functional:
                    xc_functionals.add(xc_functional)

            if property_name == 'energy_total':
                n_total_energies += 1

            if property_name == 'configuration_raw_gid':
                geometries.add(section.m_get(property_def))

            if property_name == 'section_single_configuration_calculation':
                n_calculations += 1

            if property_name == 'section_system':
                n_geometries += 1

        # Special handling for electronic/vibrational DOS and band structure:
        # these cannot currently be distinguished through the presence of a
        # single metainfo.
        searchable_quantities.discard("electronic_dos")
        searchable_quantities.discard("electronic_band_structure")
        searchable_quantities.discard("phonon_dos")
        searchable_quantities.discard("phonon_band_structure")
        if self.band_structure_electronic(entry_archive):
            searchable_quantities.add("electronic_band_structure")
        if self.band_structure_phonon(entry_archive):
            searchable_quantities.add("phonon_band_structure")
        if self.dos_electronic(entry_archive):
            searchable_quantities.add("electronic_dos")
        if self.dos_phonon(entry_archive):
            searchable_quantities.add("phonon_dos")

        self.xc_functional_names = sorted(xc_functionals)
        if len(self.xc_functional_names) > 0:
            self.xc_functional = map_functional_name_to_xc_treatment(
                get_value(self.xc_functional_names[0]))
        else:
            self.xc_functional = config.services.unavailable_value

        self.quantities = list(quantities)
        self.geometries = list(geometries)
        self.searchable_quantities = list(searchable_quantities)
        self.n_quantities = n_quantities
        self.n_calculations = n_calculations
        self.n_total_energies = n_total_energies
        self.n_geometries = n_geometries

        # grouping
        self.update_group_hash()

        # labels
        compounds = set()
        classifications = set()
        if section_system:
            for section in section_system.section_springer_material:
                compounds.update(section.springer_compound_class)
                classifications.update(section.springer_classification)

        for compound in compounds:
            self.labels.append(Label(label=compound, type='compound_class', source='springer'))
        for classification in classifications:
            self.labels.append(Label(label=classification, type='classification', source='springer'))
        self.labels_springer_compound_class = list(compounds)
        self.labels_springer_classification = list(classifications)

        aflow_id, aflow_label = None, None
        section_prototype = section_system.section_prototype if section_system else []
        if section_prototype:
            aflow_id = get_value(section_prototype[0].prototype_aflow_id)
            aflow_label = get_value(section_prototype[0].prototype_label)

        if aflow_id is not None and aflow_label is not None:
            self.labels.append(Label(label=aflow_label, type='prototype', source='aflow_prototype_library'))
            self.labels.append(Label(label=aflow_id, type='prototype_id', source='aflow_prototype_library'))

        if entry_archive.section_workflow:
            self.workflow = entry_archive.section_workflow

    def band_structure_electronic(self, entry_archive):
        """Returns whether a valid electronic band structure can be found. In
        the case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of band_k_points.
          - There is a non-empty array of band_energies.
          - The reported band_structure_kind is not "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_k_band"]
        valid = False
        for bs in self.traverse_reversed(entry_archive, path):
            kind = bs.band_structure_kind
            if kind == "vibrational" or not bs.section_k_band_segment:
                continue
            valid = True
            for segment in bs.section_k_band_segment:
                energies = segment.band_energies
                k_points = segment.band_k_points
                if not valid_array(energies) or not valid_array(k_points):
                    valid = False
                    break
            if valid:
                break
        return valid

    def dos_electronic(self, entry_archive):
        """Returns whether a valid electronic DOS can be found. In the case of
        multiple valid DOSes, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of dos_values_normalized.
          - There is a non-empty array of dos_energies.
          - The reported dos_kind is not "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_dos"]
        for dos in self.traverse_reversed(entry_archive, path):
            kind = dos.dos_kind
            energies = dos.dos_energies
            values = dos.dos_values_normalized
            if kind != "vibrational" and valid_array(energies) and valid_array(values):
                return True

        return False

    def band_structure_phonon(self, entry_archive):
        """Returns whether a valid phonon band structure can be found. In the
        case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of band_k_points.
          - There is a non-empty array of band_energies.
          - The reported band_structure_kind is "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_k_band"]
        valid = False
        for bs in self.traverse_reversed(entry_archive, path):
            kind = bs.band_structure_kind
            if kind != "vibrational" or not bs.section_k_band_segment:
                continue
            valid = True
            for segment in bs.section_k_band_segment:
                energies = segment.band_energies
                k_points = segment.band_k_points
                if not valid_array(energies) or not valid_array(k_points):
                    valid = False
                    break
            if valid:
                break

        return valid

    def dos_phonon(self, entry_archive):
        """Returns whether a valid phonon dos can be found. In the case of
        multiple valid data sources, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of dos_values_normalized.
          - There is a non-empty array of dos_energies.
          - The reported dos_kind is "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_dos"]
        for dos in self.traverse_reversed(entry_archive, path):
            kind = dos.dos_kind
            energies = dos.dos_energies
            values = dos.dos_values
            if kind == "vibrational" and valid_array(energies) and valid_array(values):
                return True

        return False

    def traverse_reversed(self, entry_archive, path):
        """Traverses the given metainfo path in reverse order. Useful in
        finding the latest reported section or value.
        """
        def traverse(root, path, i):
            sections = getattr(root, path[i])
            if isinstance(sections, list):
                for section in reversed(sections):
                    if i == len(path) - 1:
                        yield section
                    else:
                        for s in traverse(section, path, i + 1):
                            yield s
            else:
                if i == len(path) - 1:
                    yield sections
                else:
                    for s in traverse(sections, path, i + 1):
                        yield s
        for t in traverse(entry_archive, path, 0):
            if t is not None:
                yield t
