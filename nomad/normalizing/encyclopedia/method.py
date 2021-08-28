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

from typing import List
from abc import abstractmethod
from collections import OrderedDict
import numpy as np
from nomad.units import ureg

from nomad.datamodel.encyclopedia import (
    Material,
    Method,
)
from nomad.metainfo import Section
from nomad.normalizing.encyclopedia.basisset import get_basis_set
from nomad.normalizing.encyclopedia.context import Context
from nomad.utils import RestrictedDict
from nomad import config


class MethodNormalizer():
    """A base class that is used for processing method related information
    in the Encylopedia.
    """
    def __init__(self, entry_archive, logger):
        self.logger = logger
        self.entry_archive = entry_archive
        self.section_run = entry_archive.section_run[0]

    def method_id(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section):
        method_dict = RestrictedDict(
            mandatory_keys=[
                "program_name",
                "subsettings",
            ],
            forbidden_values=[None]
        )
        method_dict['program_name'] = self.section_run.program_name

        # The subclasses may define their own method properties that are to be
        # included here.
        subsettings = self.method_id_dict(method, settings_basis_set, repr_method)
        method_dict["subsettings"] = subsettings

        # If all required information is present, safe the hash
        try:
            method_dict.check(recursive=True)
        except (KeyError, ValueError) as e:
            self.logger.info("Could not create method hash: {}".format(e))
        else:
            method.method_id = method_dict.hash()

    @abstractmethod
    def method_id_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section) -> RestrictedDict:
        pass

    def group_eos_id(self, method: Method, material: Material, repr_method: Section):
        eos_dict = RestrictedDict(
            mandatory_keys=[
                "upload_id",
                "method_id",
                "formula",
            ],
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        eos_dict['upload_id'] = self.entry_archive.section_metadata.upload_id

        # Method
        eos_dict["method_id"] = method.method_id

        # The formula should be same for EoS (maybe even symmetries)
        eos_dict["formula"] = material.formula

        # Form a hash from the dictionary
        try:
            eos_dict.check(recursive=True)
        except (KeyError, ValueError) as e:
            self.logger.info("Could not create EOS hash: {}".format(e))
        else:
            method.group_eos_id = eos_dict.hash()

    def group_parametervariation_id(self, method: Method, settings_basis_set: RestrictedDict, repr_system: Section, repr_method: Section):
        # Create ordered dictionary with the values. Order is important for
        param_dict = RestrictedDict(
            mandatory_keys=[
                "upload_id",
                "program_name",
                "program_version",
                "settings_geometry",
                "subsettings",
            ],
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        param_dict['upload_id'] = self.entry_archive.section_metadata.upload_id

        # The same code and functional type is required
        param_dict['program_name'] = self.section_run.program_name
        param_dict['program_version'] = self.section_run.program_version

        # Get a string representation of the geometry. It is included as the
        # geometry should remain the same during parameter variation. By simply
        # using the atom labels and positions we assume that their
        # order/translation/rotation does not change.
        geom_dict: OrderedDict = OrderedDict()
        sec_sys = repr_system
        atom_labels = sec_sys['atom_labels']
        geom_dict['atom_labels'] = ', '.join(atom_labels)
        atom_positions = sec_sys['atom_positions']
        geom_dict['atom_positions'] = np.array2string(
            atom_positions.to(ureg.angstrom).magnitude,  # convert to Angstrom
            formatter={'float_kind': lambda x: "%.6f" % x},  # type: ignore
        ).replace('\n', '')
        cell = sec_sys['lattice_vectors']
        geom_dict['simulation_cell'] = np.array2string(
            cell.to(ureg.angstrom).magnitude,  # convert to Angstrom
            formatter={'float_kind': lambda x: "%.6f" % x},  # type: ignore
        ).replace('\n', '')
        param_dict['settings_geometry'] = geom_dict

        # The subclasses may define their own method properties that are to be
        # included here.
        subsettings = self.group_parametervariation_id_dict(method, settings_basis_set, repr_method)
        param_dict["subsettings"] = subsettings

        # Form a hash from the dictionary
        try:
            param_dict.check(recursive=True)
        except (KeyError, ValueError) as e:
            self.logger.info("Could not create parameter variation hash: {}".format(e))
        else:
            method.group_parametervariation_id = param_dict.hash()

    @abstractmethod
    def group_parametervariation_id_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section) -> RestrictedDict:
        pass

    def group_e_min(self) -> None:
        pass

    def group_type(self) -> None:
        pass

    @abstractmethod
    def normalize(self, context: Context) -> None:
        pass


class MethodDFTNormalizer(MethodNormalizer):
    """A base class that is used for processing method related information
    in the Encylopedia.
    """
    def core_electron_treatment(self, method: Method) -> None:
        treatment = config.services.unavailable_value
        code_name = self.section_run.program_name
        if code_name is not None:
            core_electron_treatments = {
                'VASP': 'pseudopotential',
                'FHI-aims': 'full all electron',
                'exciting': 'full all electron',
                'quantum espresso': 'pseudopotential'
            }
            treatment = core_electron_treatments.get(code_name, config.services.unavailable_value)
        method.core_electron_treatment = treatment

    def functional_long_name(self, method: Method, repr_method: Section) -> None:
        """'Long' form of exchange-correlation functional, list of components
        and parameters as a string: see
        https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
        """
        xc_functional = MethodDFTNormalizer.functional_long_name_from_method(repr_method, self.section_run.section_method)
        if xc_functional is config.services.unavailable_value:
            self.logger.warning(
                "Metainfo for 'XC_functional' not found, and could not "
                "compose name from 'section_XC_functionals'."
            )
        method.functional_long_name = xc_functional

    @staticmethod
    def functional_long_name_from_method(repr_method: Section, methods: List[Section]):
        """'Long' form of exchange-correlation functional, list of components
        and parameters as a string: see
        https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
        """
        linked_methods = [repr_method]
        try:
            refs = repr_method.section_method_to_method_refs
        except KeyError:
            pass
        else:
            for ref in refs:
                method_to_method_kind = ref.method_to_method_kind
                referenced_method = ref.method_to_method_ref
                if method_to_method_kind == "core_settings":
                    linked_methods.append(referenced_method)

        xc_functional = config.services.unavailable_value
        for method in linked_methods:
            try:
                section_xc_functionals = method.section_XC_functionals
            except KeyError:
                pass
            else:
                components = {}
                for component in section_xc_functionals:
                    try:
                        cname = component.XC_functional_name
                    except KeyError:
                        pass
                    else:
                        this_component = ''
                        if component.XC_functional_weight is not None:
                            this_component = str(component.XC_functional_weight) + '*'
                        this_component += cname
                        components[cname] = this_component
                result_array = []
                for name in sorted(components):
                    result_array.append(components[name])
                if len(result_array) >= 1:
                    xc_functional = '+'.join(result_array)

        return xc_functional

    def functional_type(self, method: Method) -> None:
        long_name = method.functional_long_name
        if long_name is not None:
            short_name = self.create_xc_functional_shortname(long_name)
            method.functional_type = short_name

    def method_id_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section) -> RestrictedDict:
        # Extend by DFT settings.
        hash_dict = RestrictedDict(
            mandatory_keys=(
                "functional_long_name",
                "settings_basis_set",
                "scf_threshold_energy_change",
            ),
            optional_keys=(
                "smearing_kind",
                "smearing_width",
                "number_of_eigenvalues_kpoints",
            ),
            forbidden_values=[None]
        )
        # Functional settings
        hash_dict['functional_long_name'] = method.functional_long_name

        # Basis set settings
        hash_dict['settings_basis_set'] = settings_basis_set

        # k-point sampling settings if present. Add number of kpoints as
        # detected from eigenvalues. TODO: we would like to have info on the
        # _reducible_ k-point-mesh:
        #    - grid dimensions (e.g. [ 4, 4, 8 ])
        #    - or list of reducible k-points
        smearing_kind = repr_method.smearing_kind
        if smearing_kind is not None:
            hash_dict['smearing_kind'] = smearing_kind
        smearing_width = repr_method.smearing_width
        if smearing_width is not None:
            smearing_width = '%.4f' % (smearing_width)
            hash_dict['smearing_width'] = smearing_width
        try:
            scc = self.section_run.section_single_configuration_calculation[-1]
            eigenvalues = scc.section_eigenvalues
            kpt = eigenvalues[-1].eigenvalues_kpoints
        except (KeyError, IndexError):
            pass
        else:
            if kpt is not None:
                hash_dict['number_of_eigenvalues_kpoints'] = str(len(kpt))

        # SCF convergence settings
        conv_thr = repr_method.scf_threshold_energy_change
        if conv_thr is not None:
            conv_thr = '%.13f' % (conv_thr.to(ureg.rydberg).magnitude)
            hash_dict['scf_threshold_energy_change'] = conv_thr

        return hash_dict

    def group_parametervariation_id_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section):
        """Dictionary containing the parameters used for convergence test
        grouping
        This is the source for generating the related hash."""
        param_dict = RestrictedDict(
            mandatory_keys=(
                "functional_long_name",
                "scf_threshold_energy_change",
            ),
            optional_keys=(
                "atoms_pseudopotentials",
            ),
            forbidden_values=[None]
        )

        # TODO: Add other DFT-specific properties
        # considered variations:
        #   - smearing kind/width
        #   - k point grids
        #   - basis set parameters
        # convergence threshold should be kept constant during convtest
        param_dict['functional_long_name'] = method.functional_long_name
        conv_thr = repr_method.scf_threshold_energy_change
        if conv_thr is not None:
            conv_thr = '%.13f' % (conv_thr.to(ureg.rydberg).magnitude)
        param_dict['scf_threshold_energy_change'] = conv_thr

        # Pseudopotentials are kept constant, if applicable
        if settings_basis_set is not None:
            pseudos = settings_basis_set.get('atoms_pseudopotentials', None)
            if pseudos is not None:
                param_dict['atoms_pseudopotentials'] = pseudos

        return param_dict

    def create_xc_functional_shortname(self, xc_longname):
        """Use lookup table to transform xc functional long- into shortname.
        """
        # Loof for "special" functional names listed in table
        """Easily editable table of 'short' XC functional names"""
        xc_functional_shortname = {
            'HF_X': 'HF',
            'HYB_GGA_XC_B3LYP5': 'hybrid-GGA',
            'HYB_GGA_XC_HSE06': 'hybrid-GGA',
            'BEEF-vdW': 'vdW-DF'
        }
        shortname = xc_functional_shortname.get(xc_longname, None)

        # If not, look into other options:
        if shortname is None:
            xc_functional_starts = {
                "LDA": "LDA",
                "GGA": "GGA",
                "HYB_GGA": "hybrid-GGA",
                "MGGA": "meta-GGA",
                "HYB_MGGA": "hybrid-meta-GGA",
                "HF": "HF"
            }
            sections = xc_longname.split("+")
            # decompose long name, this could be done more consistent with the
            # composition of the long name
            funcnames = []
            for section in sections:
                funcname = section.split('*')[-1]
                for func_start in xc_functional_starts:
                    if funcname.startswith(func_start):
                        funcnames.append(func_start)
                        break
            funcnames = set(funcnames)

            # Only one functional is defined
            # (usually for correlation and exchange)
            if len(funcnames) == 1:
                shortname = xc_functional_starts[func_start]
            # Two functionals that give a hybrid-GGA functional
            elif "GGA" in funcnames and "HF" in funcnames:
                shortname = "hybrid-GGA"

        if shortname is None:
            self.logger.info(
                "Could not find a functional shortname for xc_functional {}."
                .format(xc_longname)
            )

        return shortname

    def normalize(self, context: Context) -> None:
        # Fetch resources
        repr_method = context.representative_method
        repr_system = context.representative_system
        sec_enc = self.entry_archive.section_metadata.encyclopedia
        method = sec_enc.method
        material = sec_enc.material
        settings_basis_set = get_basis_set(context, self.entry_archive, self.logger)

        # Fill metainfo
        self.core_electron_treatment(method)
        self.functional_long_name(method, repr_method)
        self.functional_type(method)
        self.method_id(method, settings_basis_set, repr_method)
        self.group_eos_id(method, material, repr_method)
        self.group_parametervariation_id(method, settings_basis_set, repr_system, repr_method)


class MethodGWNormalizer(MethodDFTNormalizer):
    """A base class that is used for processing GW calculations.
    """
    def gw_starting_point(self, method: Method, repr_method: Section) -> None:
        try:
            ref = repr_method.section_method_to_method_refs[0]
            method_to_method_kind = ref.method_to_method_kind
            start_method = ref.method_to_method_ref
        except KeyError:
            pass
        else:
            if method_to_method_kind == "starting_point":
                methods = self.section_run.section_method
                xc_functional = MethodDFTNormalizer.functional_long_name_from_method(start_method, methods)
                method.gw_starting_point = xc_functional

    def functional_type(self, method: Method) -> None:
        method.functional_type = "GW"

    def gw_type(self, method: Method, repr_method: Section) -> None:
        method.gw_type = repr_method["electronic_structure_method"]

    def normalize(self, context: Context) -> None:
        # Fetch resources
        repr_method = context.representative_method
        sec_enc = self.entry_archive.section_metadata.encyclopedia
        method = sec_enc.method

        # Fill metainfo
        self.functional_type(method)
        self.gw_type(method, context.representative_method)
        self.gw_starting_point(method, repr_method)
