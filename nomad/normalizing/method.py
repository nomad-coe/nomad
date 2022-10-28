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
from abc import ABC, abstractmethod
from collections import OrderedDict
import numpy as np
from typing import Tuple, List, Union

from nomad.units import ureg
from nomad.metainfo import Section
from nomad.utils import RestrictedDict
from nomad import config
from nomad.datamodel.results import (
    Method,
    Electronic,
    Simulation,
    DFT,
    Projection,
    GW,
    xc_treatments,
)


class MethodNormalizer():
    def __init__(self, entry_archive, repr_system, material, logger):
        self.entry_archive = entry_archive
        self.repr_system = repr_system
        self.repr_method = None
        self.material = material
        self.method_name = None
        self.run = self.entry_archive.run[0]
        self.logger = logger

    def method(self) -> Method:
        """Returns a filled Method section.

        Returns:
            Filled method section.
        """
        method = Method()
        simulation = Simulation()
        repr_method = None
        method_name = config.services.unavailable_value
        methods = self.run.method
        n_methods = len(methods)

        def get_method_name(section_method):
            method_name = config.services.unavailable_value
            if section_method.electronic and section_method.electronic.method:
                method_name = section_method.electronic.method
            else:
                if section_method.gw is not None:
                    method_name = "GW"
                elif section_method.projection is not None:
                    method_name = "Projection"
            return method_name

        # If only one method is specified, use it directly
        if n_methods == 1:
            repr_method = methods[0]
            method_name = get_method_name(repr_method)
        # If several methods have been declared, we need to find the "topmost"
        # method and report it.
        elif n_methods > 1:
            # Method referencing another as "core_method". If core method was
            # given, create new merged method containing all the information.
            for sec_method in methods:
                core_method = sec_method.core_method_ref
                if core_method is not None:
                    if sec_method.electronic:
                        electronic = core_method.electronic
                        electronic = electronic if electronic else core_method.m_create(Electronic)
                        core_method.electronic.method = sec_method.electronic.method
                    repr_method = core_method
                    method_name = get_method_name(repr_method)

            # Perturbative methods: we report the "topmost" method (=method
            # that is referencing something but is not itself being
            # referenced).
            referenced_methods = set()
            for sec_method in methods:
                starting_method_ref = sec_method.starting_method_ref
                if starting_method_ref is not None:
                    referenced_methods.add(starting_method_ref.m_path())
            if len(referenced_methods) == n_methods - 1:
                for sec_method in methods:
                    if sec_method.m_path() not in referenced_methods:
                        method_name = get_method_name(sec_method)
                        if method_name != config.services.unavailable_value:
                            repr_method = sec_method
                            break
        self.repr_method = repr_method
        self.method_name = method_name
        settings_basis_set = get_basis_set(self.entry_archive, self.repr_method, self.repr_system, self.logger)
        functional_long_name = self.functional_long_name()
        method.workflow_name = self.workflow_name()

        if method_name == "GW":
            method.method_name = "GW"
            gw = GW()
            gw.type = repr_method.gw.type
            gw.starting_point = repr_method.gw.starting_point.split()
            simulation.gw = gw
        elif method_name == "Projection":
            method.method_name = "Projection"
            projection = Projection()
            if repr_method.projection.is_maximally_localized:
                projection.localization_type = 'maximally_localized'
            else:
                projection.localization_type = 'single_shot'
            simulation.projection = projection
        elif method_name in {"DFT", "DFT+U"}:
            method.method_name = "DFT"
            dft = DFT()
            dft.basis_set_type = self.basis_set_type()
            dft.basis_set_name = self.basis_set_name()
            method.method_id = self.method_id_dft(settings_basis_set, functional_long_name)
            method.parameter_variation_id = self.parameter_variation_id_dft(settings_basis_set, functional_long_name)
            dft.core_electron_treatment = self.core_electron_treatment()
            if repr_method.electronic is not None:
                if repr_method.electronic.smearing is not None:
                    dft.smearing_kind = repr_method.electronic.smearing.kind
                    dft.smearing_width = repr_method.electronic.smearing.width
                if repr_method.electronic.n_spin_channels:
                    dft.spin_polarized = repr_method.electronic.n_spin_channels > 1
                dft.van_der_Waals_method = repr_method.electronic.van_der_waals_method
                dft.relativity_method = repr_method.electronic.relativity_method
            dft.xc_functional_names = self.xc_functional_names()
            dft.xc_functional_type = self.xc_functional_type(dft.xc_functional_names)
            if repr_method.scf is not None:
                dft.scf_threshold_energy_change = repr_method.scf.threshold_energy_change
            simulation.dft = dft

        method.equation_of_state_id = self.equation_of_state_id(method.method_id, self.material.chemical_formula_hill)
        simulation.program_name = self.run.program.name
        simulation.program_version = self.run.program.version
        method.simulation = simulation
        return method

    def workflow_name(self):
        workflow_name = None
        if self.entry_archive.workflow:
            workflow_name = set(filter(lambda x: x, map(lambda x: x.type, self.entry_archive.workflow)))
        return workflow_name

    def method_id_dft(self, settings_basis_set, functional_long_name: str):
        """Creates a method id for DFT calculations if all required data is
        present.
        """
        method_dict = RestrictedDict(
            mandatory_keys=[
                "program_name",
                "functional_long_name",
                "settings_basis_set",
                "scf_threshold_energy_change",
            ],
            optional_keys=(
                "smearing_kind",
                "smearing_width",
                "number_of_eigenvalues_kpoints",
            ),
            forbidden_values=[None]
        )
        method_dict['program_name'] = self.run.program.name

        # Functional settings
        method_dict['functional_long_name'] = functional_long_name

        # Basis set settings
        method_dict['settings_basis_set'] = settings_basis_set

        # k-point sampling settings if present. Add number of kpoints as
        # detected from eigenvalues. TODO: we would like to have info on the
        # _reducible_ k-point-mesh:
        #    - grid dimensions (e.g. [ 4, 4, 8 ])
        #    - or list of reducible k-points
        try:
            smearing_kind = self.repr_method.electronic.smearing.kind
            if smearing_kind is not None:
                method_dict['smearing_kind'] = smearing_kind
            smearing_width = self.repr_method.electronic.smearing.width
            if smearing_width is not None:
                smearing_width = '%.4f' % (smearing_width)
                method_dict['smearing_width'] = smearing_width
        except Exception:
            pass

        try:
            scc = self.run.calculation[-1]
            eigenvalues = scc.eigenvalues
            kpt = eigenvalues[-1].kpoints
        except (KeyError, IndexError):
            pass
        else:
            if kpt is not None:
                method_dict['number_of_eigenvalues_kpoints'] = str(len(kpt))

        # SCF convergence settings
        try:
            conv_thr = self.repr_method.scf.threshold_energy_change
            if conv_thr is not None:
                conv_thr = '%.13f' % (conv_thr.to(ureg.rydberg).magnitude)
                method_dict['scf_threshold_energy_change'] = conv_thr
        except Exception:
            pass

        # If all required information is present, safe the hash
        try:
            method_dict.check(recursive=True)
        except (KeyError, ValueError):
            pass
        else:
            return method_dict.hash()

    def equation_of_state_id(self, method_id: str, formula: str):
        """Creates an ID that can be used to group an equation of state
        calculation found within the same upload.
        """
        eos_dict = RestrictedDict(
            mandatory_keys=[
                "upload_id",
                "method_id",
                "formula",
            ],
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        eos_dict['upload_id'] = self.entry_archive.metadata.upload_id

        # Method
        eos_dict["method_id"] = method_id

        # The formula should be same for EoS (maybe even symmetries)
        if self.material:
            eos_dict["formula"] = self.material.chemical_formula_hill

        # Form a hash from the dictionary
        try:
            eos_dict.check(recursive=True)
        except (KeyError, ValueError):
            pass
        else:
            return eos_dict.hash()

    def parameter_variation_id_dft(self, settings_basis_set, functional_long_name: str):
        """Creates an ID that can be used to group calculations that differ
        only by the used DFT parameters within the same upload.
        """
        # Create ordered dictionary with the values. Order is important for
        param_dict = RestrictedDict(
            mandatory_keys=[
                "upload_id",
                "program_name",
                "program_version",
                "settings_geometry",
                "functional_long_name",
                "scf_threshold_energy_change",
            ],
            optional_keys=(
                "atoms_pseudopotentials",
            ),
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        param_dict['upload_id'] = self.entry_archive.metadata.upload_id

        # The same code and functional type is required
        param_dict['program_name'] = self.run.program.name
        param_dict['program_version'] = self.run.program.version

        # Get a string representation of the geometry. It is included as the
        # geometry should remain the same during parameter variation. By simply
        # using the atom labels and positions we assume that their
        # order/translation/rotation does not change.
        geom_dict: OrderedDict = OrderedDict()
        try:
            atoms = self.repr_system.atoms
            atom_labels = atoms.labels
            geom_dict['atom_labels'] = ', '.join(sorted(atom_labels))
            atom_positions = atoms['positions']
            geom_dict['atom_positions'] = np.array2string(
                atom_positions.to(ureg.angstrom).magnitude,  # convert to Angstrom
                formatter={'float_kind': lambda x: "%.6f" % x},  # type: ignore
            ).replace('\n', '')
            cell = atoms['lattice_vectors']
            geom_dict['simulation_cell'] = np.array2string(
                cell.to(ureg.angstrom).magnitude,  # convert to Angstrom
                formatter={'float_kind': lambda x: "%.6f" % x},  # type: ignore
            ).replace('\n', '')
        except Exception:
            pass
        param_dict['settings_geometry'] = geom_dict

        # TODO: Add other DFT-specific properties
        # considered variations:
        #   - smearing kind/width
        #   - k point grids
        #   - basis set parameters
        # convergence threshold should be kept constant during convtest
        param_dict['functional_long_name'] = functional_long_name
        conv_thr = self.repr_method.scf.threshold_energy_change if self.repr_method.scf is not None else None
        if conv_thr is not None:
            conv_thr = '%.13f' % (conv_thr.to(ureg.rydberg).magnitude)
        param_dict['scf_threshold_energy_change'] = conv_thr

        # Pseudopotentials are kept constant, if applicable
        if settings_basis_set is not None:
            pseudos = settings_basis_set.get('atoms_pseudopotentials', None)
            if pseudos is not None:
                param_dict['atoms_pseudopotentials'] = pseudos

        # Form a hash from the dictionary
        try:
            param_dict.check(recursive=True)
        except (KeyError, ValueError):
            pass
        else:
            return param_dict.hash()

    def basis_set_type(self) -> str:
        try:
            name = self.repr_method.basis_set[0].type
        except Exception:
            name = None
        if name:
            key = name.replace('_', '').replace('-', '').replace(' ', '').lower()
            name_mapping = {
                'gaussians': 'gaussians',
                'realspacegrid': 'real-space grid',
                'planewaves': 'plane waves'
            }
            name = name_mapping.get(key, name)
        return name

    def basis_set_name(self) -> Union[str, None]:
        try:
            name = self.repr_method.basis_set[0].name
        except Exception:
            name = None
        return name

    def core_electron_treatment(self) -> str:
        treatment = config.services.unavailable_value
        code_name = self.run.program.name
        if code_name is not None:
            core_electron_treatments = {
                'VASP': 'pseudopotential',
                'FHI-aims': 'full all electron',
                'exciting': 'full all electron',
                'quantum espresso': 'pseudopotential'
            }
            treatment = core_electron_treatments.get(code_name, config.services.unavailable_value)
        return treatment

    def xc_functional_names(self) -> Union[List[str], None]:
        if self.repr_method:
            functionals = set()
            try:
                for functional_type in ['exchange', 'correlation', 'hybrid', 'contributions']:
                    functionals.update([f.name for f in self.repr_method.dft.xc_functional[functional_type]])
            except Exception:
                pass
            if functionals:
                return sorted(functionals)
        return None

    def xc_functional_type(self, xc_functionals) -> str:
        if xc_functionals:
            name = xc_functionals[0]
            return xc_treatments.get(name[:3].lower(), config.services.unavailable_value)
        else:
            return config.services.unavailable_value

    def functional_long_name(self) -> str:
        """'Long' form of exchange-correlation functional, list of components
        and parameters as a string: see
        https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
        """
        xc_functional = MethodNormalizer.functional_long_name_from_method(self.repr_method, self.run.method)
        if xc_functional is config.services.unavailable_value:
            self.logger.warning(
                "Metainfo for 'XC_functional' not found, and could not "
                "compose name from 'section_XC_functionals'."
            )
        return xc_functional

    @staticmethod
    def functional_long_name_from_method(repr_method: Section, methods: List[Section]):
        """'Long' form of exchange-correlation functional, list of components
        and parameters as a string: see
        https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
        """
        if not repr_method or not methods:
            return None

        linked_methods = [repr_method]
        if repr_method.core_method_ref is not None:
            linked_methods.append(repr_method.core_method_ref)

        xc_functional = config.services.unavailable_value
        for method in linked_methods:
            if method.dft is None:
                continue
            try:
                sec_xc_functionals = method.dft.xc_functional
            except KeyError:
                pass
            else:
                if sec_xc_functionals is None:
                    return xc_functional
                components = {}
                for functional in ['exchange', 'correlation', 'hybrid', 'contributions']:
                    for component in sec_xc_functionals[functional]:
                        try:
                            cname = component.name
                        except KeyError:
                            pass
                        else:
                            this_component = ''
                            if component.weight is not None:
                                this_component = str(component.weight) + '*'
                            this_component += cname
                            components[cname] = this_component
                    result_array = []
                    for name in sorted(components):
                        result_array.append(components[name])
                    if len(result_array) >= 1:
                        xc_functional = '+'.join(result_array)
                if method.dft.xc_functional.name is None:
                    method.dft.xc_functional.name = xc_functional

        return xc_functional


def get_basis_set(entry_archive, repr_method, repr_system, logger) -> RestrictedDict:
    """Decide which type of basis set settings are applicable to the entry and
    return a corresponding settings as a RestrictedDict.

    Args:
        entry_archive: EntryArchive from which values are extracted.
        logger: Shared logger.

    Returns:
        RestrictedDict or None: Returns the extracted settings as a
        RestrictedDict. If no suitable basis set settings could be identified,
        returns None.
    """
    settings: BasisSet = None
    program_name = entry_archive.run[0].program.name
    if program_name == "exciting":
        settings = BasisSetExciting(entry_archive, repr_method, repr_system, logger)
    elif program_name == "FHI-aims":
        settings = BasisSetFHIAims(entry_archive, repr_method, repr_system, logger)
    else:
        return None

    return settings.to_dict()


class BasisSet(ABC):
    """Abstract base class for basis set settings. The idea is to create
    subclasses that inherit this class and hierarchically add new mandatory and
    optional settings with the setup()-function.
    """
    def __init__(self, entry_archive, repr_method, repr_system, logger):
        """
        """
        self._entry_archive = entry_archive
        self._repr_method = repr_method
        self._repr_system = repr_system
        self._logger = logger
        mandatory, optional = self.setup()
        self.settings = RestrictedDict(mandatory, optional, forbidden_values=[None])

    @abstractmethod
    def to_dict(self) -> RestrictedDict:
        """Used to extract basis set settings from the archive and returning
        them as a RestrictedDict.
        """
        pass

    @abstractmethod
    def setup(self) -> Tuple:
        """Used to define a list of mandatory and optional settings for a
        subclass.

        Returns:
            Should return a tuple of two lists: the first one defining
            mandatory keys and the second one defining optional keys.
        """
        mandatory: List = []
        optional: List = []
        return mandatory, optional


class BasisSetFHIAims(BasisSet):
    """Basis set settings for 'FHI-Aims' (code-dependent).
    """
    def setup(self) -> Tuple:
        # Get previously defined values from superclass
        mandatory, optional = super().setup()

        # Add new values
        mandatory += ["fhiaims_basis"]

        return mandatory, optional

    def to_dict(self):
        # Get basis set settings for each species
        aims_bs = self._repr_method.x_fhi_aims_section_controlIn_basis_set
        if not aims_bs:
            try:
                aims_bs = self._repr_method.method_ref.x_fhi_aims_section_controlIn_basis_set
            except Exception:
                pass
        if aims_bs is not None:
            bs_by_species = {}
            for this_aims_bs in aims_bs:
                this_bs_dict = self._values_to_dict(this_aims_bs, level=2)
                this_species = this_aims_bs['x_fhi_aims_controlIn_species_name'][0]
                bs_by_species[this_species] = this_bs_dict

            # Sort alphabetically by species label
            if bs_by_species:
                basis = OrderedDict()
                for k in sorted(bs_by_species.keys()):
                    basis[k] = bs_by_species[k]
                self.settings["fhiaims_basis"] = basis

        return self.settings

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
        elif isinstance(data, (np.ndarray)):
            result = data.tolist()
        else:
            result = data
        return result

    @classmethod
    def _filtered_section_keys(cls, section):
        for k in section.keys():
            # skip JSON-specific keys
            if k == '_gIndex':
                continue
            if k == '_name':
                continue
            else:
                # json values and subsections
                yield k


class BasisSetExciting(BasisSet):
    """Basis set settings for 'Exciting' (code-dependent).
    """
    def setup(self) -> Tuple:
        # Get previously defined values from superclass
        mandatory, optional = super().setup()

        # Add new values
        mandatory += [
            "muffin_tin_settings",
            "rgkmax",
            "gkmax",
            "lo",
            "lmaxapw",
        ]

        return mandatory, optional

    def to_dict(self):
        """Special case of basis set settings for Exciting code. See list at:
        https://gitlab.mpcdf.mpg.de/nomad-lab/encyclopedia-general/wikis/FHI-visit-preparation
        """
        # Add the muffin-tin settings for each species ordered alphabetically by atom label
        try:
            groups = self._repr_system.x_exciting_section_atoms_group
            groups = sorted(groups, key=lambda group: group.x_exciting_geometry_atom_labels)
            muffin_tin_settings = OrderedDict()
            for group in groups:
                label = group.x_exciting_geometry_atom_labels
                try:
                    muffin_tin_settings["{}_muffin_tin_radius".format(label)] = "%.6f" % (group.x_exciting_muffin_tin_radius.to(ureg.angstrom).magnitude)
                except Exception:
                    muffin_tin_settings["{}_muffin_tin_radius".format(label)] = None
                try:
                    muffin_tin_settings["{}_muffin_tin_points".format(label)] = "%d" % group.x_exciting_muffin_tin_points
                except Exception:
                    muffin_tin_settings["{}_muffin_tin_points".format(label)] = None
            self.settings["muffin_tin_settings"] = muffin_tin_settings
        except Exception:
            pass

        # Other important method settings
        system = self._repr_system
        try:
            self.settings['rgkmax'] = "%.6f" % (system.x_exciting_rgkmax.magnitude)
        except Exception:
            pass
        try:
            self.settings['gkmax'] = "%.6f" % (1e-10 * system.x_exciting_gkmax.magnitude)
        except Exception:
            pass
        try:
            self.settings['lo'] = "%d" % (system.x_exciting_lo)
        except Exception:
            pass
        try:
            self.settings['lmaxapw'] = "%d" % (system.x_exciting_lmaxapw)
        except Exception:
            pass

        return self.settings
