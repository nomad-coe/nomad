# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import OrderedDict
from abc import ABC, abstractmethod
from nomad.parsing.backend import Section
import re
import numpy

J_to_Ry = 4.587425e+17


class SettingsBasisSet(ABC):
    """Abstract base class for basis set settings

    Provides a factory() static method for delegating to the concrete
    implementation applicable for a given calculation.
    """
    def __init__(self, context, backend, logger):
        """
        """
        self.context = context
        self.backend = backend
        self.logger = logger

    @staticmethod
    def factory(ctx, backend, logger):
        """Decide which type of basis set settings are applicable to the entry
        and return a corresponding SettingsBasisSet class.
        """
        # Check if there is a code-specific SettingsBasisSet
        for cls in SettingsBasisSetCodeDependent.__subclasses__():
            if cls.is_basis_set_for(backend):
                return cls(ctx, backend, logger)

        # Check if there is a generic SettingsBasisSet
        for cls in SettingsBasisSetGeneric.__subclasses__():
            if cls.is_basis_set_for(backend):  # pylint: disable=E1120
                return cls(ctx, backend, logger)

        # Raise Exception in case we did not find an implementation
        raise ValueError("No SettingsBasisSet available for this entry.")

    @staticmethod
    @abstractmethod
    def is_basis_set_for(backend):
        """Report if this (sub-)class is applicable to the given calculation"""
        return False

    def to_dict(self, result_dict=None):
        """Dictionary representation of basis set settings; optionally, add
        as new entries to a pre-existing dictionary 'result_dict'"""
        # create new dictionary if necessary
        if result_dict is None:
            result_dict = OrderedDict()
        return result_dict


class SettingsBasisSetGeneric(SettingsBasisSet):
    """Base class for code-independent / generic basis set settings"""
    pass


class SettingsBasisSetCodeDependent(SettingsBasisSet):
    """Base class for code-dependent basis set settings"""
    pass


class SettingsBasisSetAuxiliary(SettingsBasisSet):
    """Base class for auxiliary data; to be inherited by other classes
    """
    def is_basis_set_for(self, backend):
        # auxiliary data should never be autodetected as basis set
        return False


class SettingsBasisSetAuxiliary_Pseudopotentials(SettingsBasisSetAuxiliary):
    """Utility class for any calculations employing pseudopotentials
    """
    def to_dict(self, result_dict=None):
        """Extract the name of the pseudopotential used for each atom,
        and return the resulting list in dictionary key
        'atoms_pseudopotentials'
        """
        result_dict = super().to_dict(result_dict)
        atom_kind_list = self.backend.get('section_method_atom_kind')
        if atom_kind_list is None:
            return result_dict

        # Default lookup strategy is by atom label
        #   build lookup table: label -> atom_kind_index
        atom_kind_by_label = {}
        for atom_kind_i in range(len(atom_kind_list)):
            labels = atom_kind_list[atom_kind_i].get('method_atom_kind_label')
            if labels is not None:
                for label in labels:
                    atom_kind_by_label[label] = atom_kind_list[atom_kind_i]
        sec_sys = self.context.representative_system
        atom_labels = sec_sys.get('atom_labels')
        atom_kinds = [None] * len(atom_labels)
        for atom_i in range(len(atom_labels)):
            atom_kinds[atom_i] = atom_kind_by_label.get(
                atom_labels[atom_i], None)

        # Fallback for VASP, which provides atom->index ref list but no
        #   'method_atom_kind_label'
        sec_method = self.context.representative_method
        atom_kind_index = sec_method.get('x_vasp_atom_kind_refs')
        if atom_kind_index is not None:
            for atom_i in range(len(atom_kind_index)):
                kind_i = atom_kind_index[atom_i]
                if kind_i < len(atom_kind_list):
                    atom_kinds[atom_i] = atom_kind_list[kind_i]
        result = []
        have_names = 0
        for atom_kind in atom_kinds:
            if atom_kind is None:
                result.append(None)
                continue
            pp_name = atom_kind.get('method_atom_kind_pseudopotential_name')
            if pp_name is not None:
                pp_name = pp_name[0]
            elif self.backend['program_name'] == "Quantum Espresso":
                # Fallback for quantum espresso code-specific metaInfo
                fname = atom_kind.get('x_qe_pp_filename')
                if fname is not None:
                    # Use basename of PP file
                    m = re.match(r".*?([^/]+)$", fname[0])
                    pp_name = m.group(1)
            if pp_name is not None:
                have_names += 1
            result.append(pp_name)
        if have_names > 0:
            result_dict['atoms_pseudopotentials'] = result
        return result_dict


class SettingsBasisSetGenericPlaneWaves(
    SettingsBasisSetGeneric,
    SettingsBasisSetAuxiliary_Pseudopotentials,
):
    """Basis set settings for plane-waves codes
    """
    @staticmethod
    def is_basis_set_for(backend):
        return backend.get('program_basis_set_type') == 'plane waves'

    def to_dict(self, result_dict=None):
        result_dict = super().to_dict(result_dict)
        try:
            result_dict['plane_wave_cutoff_wavefunction'] = '%.4f' % (
                self._plane_wave_cutoff('wavefunction'))
        except Exception:
            result_dict['plane_wave_cutoff_wavefunction'] = None

        try:
            result_dict['plane_wave_cutoff_density'] = '%.4f' % (
                self._plane_wave_cutoff('density'))
        except Exception:
            # NOTE: treat charge-density cutoff for now as optional
            pass

        return result_dict

    def _plane_wave_cutoff(self, kind):
        """get plane-wave cutoff for a kind of grid; presently defined kinds
        are 'wavefunction' and 'density'.
        Returns value in Rydberg units."""
        meth_bs = self.backend.get('section_method_basis_set')
        if meth_bs is None:
            return None
        for component in meth_bs:
            component_kind = component.get('method_basis_set_kind')
            if component_kind != kind:
                continue
            bs_id = component.get(
                'mapping_section_method_basis_set_cell_associated')
            bs_ref = 'section_basis_set_cell_dependent:%d' % (bs_id)
            bs = self.backend.get(bs_ref)
            cutoff = bs.get('basis_set_planewave_cutoff')[0] * J_to_Ry
            return cutoff


class SettingsBasisSetCodeDependentFhiAims(
    SettingsBasisSetCodeDependent,
):
    """Basis set settings for 'FHI-Aims' (code-dependent)
    """
    @staticmethod
    def is_basis_set_for(backend):
        return backend.get('program_name') == 'FHI-aims'

    def to_dict(self, result_dict=None):
        """Special case of basis set settings for FHI-Aims code.
        """
        result_dict = super().to_dict(result_dict)
        result_dict['FhiAims_basis'] = None

        # get section_method/x_fhi_aims_section_controlIn_basis_set
        #   lookup via shortcut setup in Nomad/Archive/calculation.py
        aims_bs = self.backend.get('x_fhi_aims_section_controlIn_basis_set')
        if aims_bs is None:
            self.logger.warning(
                "could not resolve x_fhi_aims_section_controlIn_basis_set"
            )
            return result_dict

        # each 'x_fhi_aims_section_controlIn_basis_set' describes basis for
        # a species; store per-species so we can sort later on
        bs_by_species = {}
        for this_aims_bs in aims_bs:
            this_bs_dict = self._values_to_dict(this_aims_bs, level=2)
            this_species = this_aims_bs['x_fhi_aims_controlIn_species_name'][0]
            if this_species in bs_by_species:
                self.logger.warning(
                    "multiple definitions of x_fhi_aims_section_controlIn_basis_set for species " + this_species
                )
                return result_dict
            bs_by_species[this_species] = this_bs_dict

        # result: sorted alphabetically by species label
        if bs_by_species:
            result_dict['FhiAims_basis'] = OrderedDict()
            for k in sorted(bs_by_species.keys()):
                result_dict['FhiAims_basis'][k] = bs_by_species[k]
        return result_dict

    # TODO: this should go as .keys() method to archive.py hdf5/json
    @classmethod
    def _filtered_section_keys(cls, section):
        for k in section.keys():
            # skip json-specific keys
            if k == 'gIndex':
                continue
            if k == 'name':
                continue
            if k == 'references':
                continue
            if k == 'type':
                continue
            # skip hdf5-specific keys
            if k.endswith("-index"):
                continue
            if k.endswith("-v"):
                # hdf5 values
                yield k[:-2]
            else:
                # json values and subsections
                yield k

    @classmethod
    def _values_to_dict(cls, data, level=0):
        result = None
        if data is None:
            return None
        elif isinstance(data, (Section, dict)):
            result = OrderedDict()
            for k in sorted(cls._filtered_section_keys(data)):
                v = data.get(k, None)
                result[k] = cls._values_to_dict(v, level=level + 1)
        elif isinstance(data, (list)):
            result = []
            for k in range(len(data)):
                v = data[k]
                result.append(cls._values_to_dict(v, level=level + 1))
        elif isinstance(data, (numpy.ndarray)):
            result = data.tolist()
        else:
            result = data
        return result


class SettingsBasisSetCodeDependent_Exciting(
    SettingsBasisSetCodeDependent,
):
    """Basis set settings for 'Exciting' (code-dependent)
    """
    @staticmethod
    def is_basis_set_for(backend):
        return backend.get('program_name') == 'exciting'

    def to_dict(self, result_dict=None):
        """Special case of basis set settings for Exciting code. See list at:
        https://gitlab.mpcdf.mpg.de/nomad-lab/encyclopedia-general/wikis/FHI-visit-preparation
        """
        result_dict = super().to_dict(result_dict)

        system = self.context.representative_system

        try:
            result_dict['muffin_tin_radius'] = ', '.join(map(
                lambda r: "%.6f" % (r),
                system['x_exciting_muffin_tin_radius'] * 1e+10))
        except Exception:
            result_dict['muffin_tin_radius'] = None

        try:
            result_dict['rgkmax'] = "%.6f" % (
                system['x_exciting_rgkmax'][0] * 1e+10)
        except Exception:
            result_dict['rgkmax'] = None

        try:
            result_dict['gkmax'] = "%.6f" % (
                system['x_exciting_gkmax'][0] * 1e-10)
        except Exception:
            result_dict['gkmax'] = None

        try:
            result_dict['lo'] = "%d" % (system['x_exciting_lo'][0])
        except Exception:
            result_dict['lo'] = None

        try:
            result_dict['lmaxapw'] = "%d" % (system['x_exciting_lmaxapw'][0])
        except Exception:
            result_dict['lmaxapw'] = None

        try:
            result_dict['muffin_tin_points'] = ', '.join(map(
                lambda r: "%d" % (r),
                system['x_exciting_muffin_tin_points']))
        except Exception:
            result_dict['muffin_tin_points'] = None

        return result_dict
