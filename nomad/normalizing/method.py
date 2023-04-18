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
from ase.dft.kpoints import monkhorst_pack, get_monkhorst_pack_size_and_offset
from collections import OrderedDict
import re
import numpy as np
from typing import List, Tuple, Union

from nomad.units import ureg
from nomad.metainfo import Section
from nomad.metainfo.util import MTypes
from nomad.utils import RestrictedDict
from nomad import config
from nomad.datamodel.metainfo.simulation.method import KMesh
from nomad.datamodel.results import (
    Method,
    Electronic,
    Simulation,
    HubbardKanamoriModel,
    DFT,
    Projection,
    GW,
    DMFT,
    Precision,
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
        '''Returns the populated results.method section.'''
        method = Method()
        simulation = Simulation()
        repr_method = None
        method_name = config.services.unavailable_value
        workflow_name = None
        methods = self.run.method
        n_methods = len(methods)

        def get_method_name(section_method):
            method_name = config.services.unavailable_value
            if section_method.m_xpath('electronic') and section_method.electronic.m_xpath('method'):
                method_name = section_method.electronic.method
            elif section_method.m_xpath('gw') is not None:
                method_name = 'GW'
            elif section_method.m_xpath('projection') is not None:
                method_name = 'Projection'
            elif section_method.m_xpath('dmft') is not None:
                method_name = 'DMFT'
            elif section_method.m_xpath('core_hole') is not None:
                method_name = 'CoreHole'
            elif section_method.m_xpath('bse') is not None:
                method_name = 'BSE'
            return method_name

        def functional_long_name_from_method(methods):
            '''
            Long form of exchange-correlation functional, list of components and parameters
            as a string: see https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
            '''
            if not self.repr_method or not methods:
                return None

            linked_methods = [self.repr_method]
            xc_functional = config.services.unavailable_value
            for method in linked_methods:
                sec_xc_functionals = None
                if method.dft:
                    sec_xc_functionals = method.dft.xc_functional

                if sec_xc_functionals:
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
                    if method.dft and method.dft.xc_functional.name is None:
                        method.dft.xc_functional.name = xc_functional
            return xc_functional

        def get_basis_set() -> RestrictedDict:
            '''
            Decide which type of basis set settings are applicable to the entry and
            return a corresponding settings as a RestrictedDict.

            Returns:
                RestrictedDict or None: Returns the extracted settings as a
                RestrictedDict. If no suitable basis set settings could be identified,
                returns None.
            '''
            settings: BasisSet = None
            program_name = None
            if len(self.entry_archive.run) > 0 and self.run.program:
                program_name = self.run.program.name
            if program_name == 'exciting':
                settings = BasisSetExciting(self.entry_archive, self.repr_method, self.repr_system, self.logger)
            elif program_name == 'FHI-aims':
                settings = BasisSetFHIAims(self.entry_archive, self.repr_method, self.repr_system, self.logger)
            else:
                return None

            return settings.to_dict()

        # workflow_name
        if self.entry_archive.workflow:
            workflow_name = list(filter(lambda x: x, map(lambda x: x.type, self.entry_archive.workflow)))[0]
        method.workflow_name = workflow_name

        # repr_method and method_name
        #   If workflow_name is 'GW', set repr_method as the corresponding method_ref
        if method.workflow_name and method.workflow_name == 'GW':
            try:
                dft_method = self.entry_archive.workflow[0].workflows_ref[0].calculations_ref[-1].method_ref
                gw_method = self.entry_archive.workflow[0].workflows_ref[1].calculations_ref[-1].method_ref
                repr_method = dft_method
                method_name = f'{get_method_name(repr_method)}+GW'
            except Exception:
                self.logger.warning('Error finding the DFT and GW method sections from workflow refs.')
                return method
        #   If only one method is specified, use it directly
        elif n_methods == 1:
            repr_method = methods[0]
            method_name = get_method_name(repr_method)
        #   If several methods have been declared, we need to find the "topmost" method and report it.
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
        if 'DFT' in self.method_name:
            functional_long_name = functional_long_name_from_method(self.run.method)
            settings_basis_set = get_basis_set()

        # Populating method metainfo
        if self.method_name in ['DFT', 'DFT+U']:
            method.method_name = 'DFT'
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
            try:
                dft.xc_functional_names = self.xc_functional_names(self.repr_method.dft.xc_functional)
                dft.xc_functional_type = self.xc_functional_type(dft.xc_functional_names)
                dft.exact_exchange_mixing_factor = self.exact_exchange_mixing_factor(dft.xc_functional_names)
            except Exception:
                self.logger.warning('Error extracting the DFT XC functional names.')
            if repr_method.scf is not None:
                dft.scf_threshold_energy_change = repr_method.scf.threshold_energy_change
            simulation.dft = dft
            hubbard_kanamori_models = self.hubbard_kanamori_model(methods)
            simulation.dft.hubbard_kanamori_model = hubbard_kanamori_models if len(hubbard_kanamori_models) else None
        elif method_name == 'GW':
            method.method_name = 'GW'
            gw = GW()
            gw.type = self.repr_method.gw.type
            simulation.gw = gw
        elif method_name == 'DFT+GW':
            method.method_name = 'GW'
            gw = GW()
            gw.type = gw_method.gw.type
            try:
                gw.starting_point_names = self.xc_functional_names(dft_method.dft.xc_functional)
                gw.starting_point_type = self.xc_functional_type(gw.starting_point_names)
            except Exception:
                self.logger.warning('Error extracting the GW starting point names.')
            gw.basis_set_type = self.basis_set_type()
            gw.basis_set_name = self.basis_set_name()
            simulation.gw = gw
        elif method_name == 'Projection':
            method.method_name = 'Projection'
            projection = Projection()
            if hasattr(self.repr_method.projection, 'wannier'):
                projection.type = 'wannier'
                if self.repr_method.projection.wannier.is_maximally_localized:
                    projection.localization_type = 'maximally_localized'
                else:
                    projection.localization_type = 'single_shot'
            elif hasattr(self.repr_method.projection, 'slater_koster'):
                projection.type = 'slater_koster'
            else:
                projection.type = 'custom'
            simulation.projection = projection
        elif method_name == 'DMFT':
            method.method_name = 'DMFT'
            dmft = DMFT()
            dmft.impurity_solver_type = self.repr_method.dmft.impurity_solver
            dmft.total_filling = 0.5 * np.sum(self.repr_method.dmft.n_correlated_electrons) / np.sum(repr_method.dmft.n_correlated_orbitals)
            dmft.inverse_temperature = self.repr_method.dmft.inverse_temperature
            dmft.magnetic_state = self.repr_method.dmft.magnetic_state
            # TODO update U to be U/W when linking between DFT>Projection>DMFT (@ should
            # be extracted from the bands obtained by Projection).
            model_hamiltonian = self.repr_method.starting_method_ref.lattice_model_hamiltonian
            if model_hamiltonian is not None:
                dmft.u = model_hamiltonian[0].hubbard_kanamori_model[0].u  # taking U,JH values from the first atom
                if dmft.u.magnitude != 0.0:
                    dmft.hunds_hubbard_ratio = model_hamiltonian[0].hubbard_kanamori_model[0].jh.magnitude / dmft.u.magnitude
            simulation.dmft = dmft
        elif method_name == 'CoreHole':
            method.method_name = 'CoreHole'
        elif method_name == 'BSE':
            method.method_name = 'BSE'

        # Fill meshes
        if self.run.m_xpath('method[-1].frequency_mesh'):
            freq_mesh = self.run.method[-1].frequency_mesh
            freq_mesh.dimensionality = 1 if freq_mesh.dimensionality is None else freq_mesh.dimensionality

        if self.run.m_xpath('method[-1].time_mesh'):
            time_mesh = self.run.method[-1].time_mesh
            time_mesh.dimensionality = 1 if time_mesh.dimensionality is None else time_mesh.dimensionality

        if self.run.m_xpath('method[-1].k_mesh'):
            k_mesh = self.run.method[-1].k_mesh
            k_mesh.dimensionality = 3 if not k_mesh.dimensionality else k_mesh.dimensionality
            # Normalize k mesh from grid sampling
            if k_mesh.grid is not None:
                k_mesh.n_points = np.prod(k_mesh.grid) if not k_mesh.n_points else k_mesh.n_points
                if k_mesh.sampling_method == 'Gamma-centered':
                    k_mesh.points = np.meshgrid(*[np.linspace(0, 1, n) for n in k_mesh.grid])  # this assumes a gamma-centered grid: we really need the `sampling_method` to be sure
                elif k_mesh.sampling_method == 'Monkhorst-Pack':
                    try:
                        k_mesh.points += monkhorst_pack(k_mesh.grid)
                    except ValueError:
                        pass  # this is a quick workaround: k_mesh.grid should be symmetry reduced
        else:
            if self.run.m_xpath('calculation[-1].eigenvalues[-1].kpoints') is not None:
                k_mesh = self.run.method[-1].m_create(KMesh)
                k_mesh.points = self.run.calculation[-1].eigenvalues[-1].kpoints
                k_mesh.n_points = len(k_mesh.points) if not k_mesh.n_points else k_mesh.n_points
                k_mesh.grid = [len(set(k_mesh.points[:, i])) for i in range(3)]
                if not k_mesh.sampling_method:
                    try:  # TODO doublecheck
                        _, k_grid_offset = get_monkhorst_pack_size_and_offset(k_mesh.points)
                        if not k_grid_offset.all():
                            k_mesh.sampling_method = 'Monkhorst-Pack'
                    except ValueError:
                        k_mesh.sampling_method = 'Gamma-centered'

        # Fill the presicion section
        k_lattices = self.run.m_xpath('system[-1].atoms.lattice_vectors_reciprocal')
        grid = self.run.m_xpath('method[-1].k_mesh.grid')
        k_line_density = self.calc_k_line_density(k_lattices, grid)
        if k_line_density:
            if not simulation.precision:
                simulation.precision = Precision()
            if not simulation.precision.k_line_density:
                simulation.precision.k_line_density = k_line_density

        method.equation_of_state_id = self.equation_of_state_id(method.method_id, self.material.chemical_formula_hill)
        simulation.program_name = self.run.program.name
        simulation.program_version = self.run.program.version
        method.simulation = simulation
        return method

    def calc_k_line_density(self, k_lattices: List[List[float]],
                            nks: List[int]) -> Union[float, None]:
        '''
        Compute the lowest k_line_density value:
        k_line_density (for a uniformly spaced grid) is the number of k-points per reciprocal length unit
        '''
        # Check consistency of input
        try:
            if len(k_lattices) != 3 or len(nks) != 3:
                return None
        except (TypeError, ValueError):
            return None

        # Compute k_line_density
        struc_type = self.material.structural_type
        if struc_type == 'bulk':
            return min([nk / (np.linalg.norm(k_lattice))
                        for k_lattice, nk in zip(k_lattices, nks)])
        else:
            return None

    def hubbard_kanamori_model(self, methods) -> List[HubbardKanamoriModel]:
        '''Generate a list of normalized HubbardKanamoriModel for `results.method`'''
        hubbard_kanamori_models = []
        for sec_method in methods:
            for param in sec_method.atom_parameters:
                if param.hubbard_kanamori_model is not None:
                    hubb_run = param.hubbard_kanamori_model
                    if all([hubb_run.orbital, hubb_run.double_counting_correction]):
                        hubb_results = HubbardKanamoriModel()
                        hubb_results.atom_label = param.label
                        valid = False
                        for quant in hubb_run.m_def.quantities:
                            quant_value = getattr(hubb_run, quant.name)
                            # all false values, including zero are ignored
                            if quant_value:
                                setattr(hubb_results, quant.name, quant_value)
                                if quant.type in MTypes.num:
                                    valid = True
                        # do not save if all parameters are set at 0
                        if not valid: continue
                        # U_effective technically only makes sense for Dudarev
                        # but it is computed anyhow to act as a trigger for DFT+U
                        if hubb_results.u_effective is None and hubb_results.u is not None:
                            hubb_results.u_effective = hubb_results.u
                            if hubb_results.j is not None:
                                hubb_results.u_effective -= hubb_results.j
                        hubbard_kanamori_models.append(hubb_results)
        return hubbard_kanamori_models

    def method_id_dft(self, settings_basis_set, functional_long_name: str):
        '''Creates a method id for DFT calculations if all required data is
        present.
        '''
        method_dict = RestrictedDict(
            mandatory_keys=[
                'program_name',
                'functional_long_name',
                'settings_basis_set',
                'scf_threshold_energy_change',
            ],
            optional_keys=(
                'smearing_kind',
                'smearing_width',
                'number_of_eigenvalues_kpoints',
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
        '''Creates an ID that can be used to group an equation of state
        calculation found within the same upload.
        '''
        eos_dict = RestrictedDict(
            mandatory_keys=[
                'upload_id',
                'method_id',
                'formula',
            ],
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        eos_dict['upload_id'] = self.entry_archive.metadata.upload_id

        # Method
        eos_dict['method_id'] = method_id

        # The formula should be same for EoS (maybe even symmetries)
        if self.material:
            eos_dict['formula'] = self.material.chemical_formula_hill

        # Form a hash from the dictionary
        try:
            eos_dict.check(recursive=True)
        except (KeyError, ValueError):
            pass
        else:
            return eos_dict.hash()

    def parameter_variation_id_dft(self, settings_basis_set, functional_long_name: str):
        '''Creates an ID that can be used to group calculations that differ
        only by the used DFT parameters within the same upload.
        '''
        # Create ordered dictionary with the values. Order is important for
        param_dict = RestrictedDict(
            mandatory_keys=[
                'upload_id',
                'program_name',
                'program_version',
                'settings_geometry',
                'functional_long_name',
                'scf_threshold_energy_change',
            ],
            optional_keys=(
                'atoms_pseudopotentials',
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

    def xc_functional_names(self, method_xc_functional: Section) -> Union[List[str], None]:
        if self.repr_method:
            functionals = set()
            try:
                for functional_type in ['exchange', 'correlation', 'hybrid', 'contributions']:
                    functionals.update([f.name for f in method_xc_functional[functional_type]])
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

    def exact_exchange_mixing_factor(self, xc_functional_names):
        '''Assign the exact exachange mixing factor to `results` section when explicitly stated.
        Else, fall back on XC functional default.'''
        def scan_patterns(patterns, xc_name) -> bool:
            return any(x for x in patterns if re.search('_' + x + '$', xc_name))

        if self.repr_method.dft:
            xc_functional = self.repr_method.dft.xc_functional
            for hybrid in xc_functional.hybrid:
                if hybrid.parameters:
                    if 'exact_exchange_mixing_factor' in hybrid.parameters.keys():
                        return hybrid.parameters['exact_exchange_mixing_factor']
        for xc_name in xc_functional_names:
            if not re.search('_XC?_', xc_name): continue
            if re.search('_B3LYP[35]?$', xc_name): return .2
            elif scan_patterns(['HSE', 'PBEH', 'PBE_MOL0', 'PBE_SOL0'], xc_name): return .25
            elif re.search('_M05$', xc_name): return .28
            elif re.search('_PBE0_13$', xc_name): return 1 / 3
            elif re.search('_PBE38$', xc_name): return 3 / 8
            elif re.search('_PBE50$', xc_name): return .5
            elif re.search('_M06_2X$', xc_name): return .54
            elif scan_patterns(['M05_2X', 'PBE_2X'], xc_name): return .56


class BasisSet(ABC):
    '''Abstract base class for basis set settings. The idea is to create
    subclasses that inherit this class and hierarchically add new mandatory and
    optional settings with the setup()-function.
    '''
    def __init__(self, entry_archive, repr_method, repr_system, logger):
        self._entry_archive = entry_archive
        self._repr_method = repr_method
        self._repr_system = repr_system
        self._logger = logger
        mandatory, optional = self.setup()
        self.settings = RestrictedDict(mandatory, optional, forbidden_values=[None])

    @abstractmethod
    def to_dict(self) -> RestrictedDict:
        '''Used to extract basis set settings from the archive and returning
        them as a RestrictedDict.
        '''
        pass

    @abstractmethod
    def setup(self) -> Tuple:
        '''Used to define a list of mandatory and optional settings for a
        subclass.

        Returns:
            Should return a tuple of two lists: the first one defining
            mandatory keys and the second one defining optional keys.
        '''
        mandatory: List = []
        optional: List = []
        return mandatory, optional


class BasisSetFHIAims(BasisSet):
    '''Basis set settings for FHI-Aims (code-dependent).
    '''
    def setup(self) -> Tuple:
        # Get previously defined values from superclass
        mandatory, optional = super().setup()

        # Add new values
        mandatory += ['fhiaims_basis']

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
                self.settings['fhiaims_basis'] = basis

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
    '''Basis set settings for Exciting (code-dependent).
    '''
    def setup(self) -> Tuple:
        # Get previously defined values from superclass
        mandatory, optional = super().setup()

        # Add new values
        mandatory += [
            'muffin_tin_settings',
            'rgkmax',
            'gkmax',
            'lo',
            'lmaxapw',
        ]

        return mandatory, optional

    def to_dict(self):
        '''Special case of basis set settings for Exciting code. See list at:
        https://gitlab.mpcdf.mpg.de/nomad-lab/encyclopedia-general/wikis/FHI-visit-preparation
        '''
        # Add the muffin-tin settings for each species ordered alphabetically by atom label
        try:
            groups = self._repr_system.x_exciting_section_atoms_group
            groups = sorted(groups, key=lambda group: group.x_exciting_geometry_atom_labels)
            muffin_tin_settings = OrderedDict()
            for group in groups:
                label = group.x_exciting_geometry_atom_labels
                try:
                    muffin_tin_settings['{}_muffin_tin_radius'.format(label)] = "%.6f" % (group.x_exciting_muffin_tin_radius.to(ureg.angstrom).magnitude)
                except Exception:
                    muffin_tin_settings['{}_muffin_tin_radius'.format(label)] = None
                try:
                    muffin_tin_settings['{}_muffin_tin_points'.format(label)] = "%d" % group.x_exciting_muffin_tin_points
                except Exception:
                    muffin_tin_settings['{}_muffin_tin_points'.format(label)] = None
            self.settings['muffin_tin_settings'] = muffin_tin_settings
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
