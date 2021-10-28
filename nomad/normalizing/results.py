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

import json
import re
import numpy as np
from typing import Dict, List, Union, Any
import ase.data

from nomad import config
from nomad import atomutils
from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.encyclopedia import EncyclopediaMetadata
from nomad.datamodel.optimade import OptimadeEntry, Species
from nomad.datamodel.metainfo.simulation.system import System, Symmetry as SystemSymmetry
from nomad.datamodel.metainfo.simulation.method import Electronic
from nomad.datamodel.results import (
    BandGap,
    Results,
    Material,
    Method,
    GeometryOptimizationProperties,
    GeometryOptimizationMethod,
    Properties,
    SpectroscopyProperties,
    Symmetry,
    Structures,
    Structure,
    StructureOriginal,
    StructurePrimitive,
    StructureConventional,
    StructureOptimized,
    LatticeParameters,
    WyckoffSet,
    Simulation,
    DFT,
    GW,
    EnergyVolumeCurve,
    BulkModulus,
    ShearModulus,
    MechanicalProperties,
    xc_treatments,
    ElectronicProperties,
    VibrationalProperties,
    BandStructureElectronic,
    BandStructurePhonon,
    DOSElectronic,
    DOSPhonon,
    EnergyFreeHelmholtz,
    HeatCapacityConstantVolume,
    Experiment,
    EELS,
)

re_label = re.compile("^([a-zA-Z][a-zA-Z]?)[^a-zA-Z]*")
elements = set(ase.data.chemical_symbols)


def label_to_chemical_symbol(label: str) -> Union[str, None]:
    """Tries to extract a valid chemical symbol from a label. Currently can
    only handle labels that correspond to chemical names (no matter what case)
    and that possbily have a non-alphabetic postfix. Raises an error if no
    match can be made.
    """
    match = re_label.match(label)
    symbol = None
    if match:
        test = match.group(1).capitalize()
        if test in elements:
            symbol = test
    if symbol is None:
        raise ValueError(
            "Could not identify a chemical element from the given label: {}".format(label)
        )
    return symbol


def valid_array(array: Any) -> bool:
    """Checks if the given variable is a non-empty array.
    """
    return array is not None and len(array) > 0


def isint(value: Any) -> bool:
    """Checks if the given variable can be interpreted as an integer.
    """
    try:
        int(value)
        return True
    except ValueError:
        return False


class ResultsNormalizer(Normalizer):
    domain = None

    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        results = self.entry_archive.results
        if results is None:
            results = self.entry_archive.m_create(Results)
        if results.properties is None:
            results.m_create(Properties)

        if self.section_run:
            self.normalize_run(logger=self.logger)

        if self.entry_archive.section_measurement and len(self.entry_archive.section_measurement) > 0:
            self.normalize_measurement(self.entry_archive.section_measurement[0], logger=self.logger)

        # Add the present quantities. The full path will be saved for each
        # property, and if the leaf quantity/section name is unambiguous, also
        # it will be saved. Each repeating section will be investigated as well
        # to find all present properties.
        available_properties = set()
        for section, m_def, _ in results.properties.m_traverse():
            parent_path = section.m_path()
            path = "{}/{}".format(parent_path if parent_path != "/" else "", m_def.name)[20:]
            parts = filter(lambda x: not isint(x), path.split("/"))
            path = ".".join(parts)
            available_properties.add(path)
        shorthand_prop = list()
        for prop in available_properties:
            name = prop.rsplit(".", 1)[-1]
            shorthand_prop.append(name)
        u, c = np.unique(shorthand_prop, return_counts=True)
        shorthand_prop = u[c == 1]
        available_properties |= set([str(x) for x in shorthand_prop])
        results.properties.available_properties = list(available_properties)

    def normalize_measurement(self, measurement, logger) -> None:
        results = self.entry_archive.results
        experiment = measurement.section_metadata.section_experiment

        # Method
        if results.method is None:
            results.m_create(Method)
        if results.method.experiment is None:
            results.method.m_create(Experiment)
        method_name = experiment.method_name
        if method_name == 'electron energy loss spectroscopy':
            try:
                device_settings = measurement.section_metadata.section_instrument.section_device_settings
                n_settings = len(device_settings)
                if n_settings > 1:
                    raise ValueError("Found multiple device settings.")
                device_settings = device_settings[0]
            except Exception as e:
                logger.error('no unique device settings available', exc_info=e)
            else:
                results.method.method_name = 'EELS'
                eels = results.method.experiment.m_create(EELS)
                eels.resolution = device_settings.resolution
                eels.detector_type = device_settings.detector_type
                eels.min_energy = device_settings.min_energy
                eels.max_energy = device_settings.max_energy
        elif method_name == 'XPS':
            results.method.method_name = 'XPS'
        else:
            logger.error('unknown measurement method', data=results.method.method_name)

        # Material
        if results.material is None:
            results.m_create(Material)
        try:
            material = measurement.section_metadata.section_sample.section_material[0]
            results.material.elements = material.elements if material.elements else []
            atoms = None
            if material.formula:
                results.material.chemical_formula_descriptive = material.formula
                try:
                    atoms = ase.Atoms(material.formula)
                except Exception as e:
                    logger.warn('could not normalize formula, using elements next', exc_info=e)

            if atoms is None:
                atoms = ase.Atoms(''.join(material.elements))

            results.material.chemical_formula_descriptive = atoms.get_chemical_formula(mode='hill')
            results.material.chemical_formula_reduced = atoms.get_chemical_formula(mode='reduce')
            results.material.chemical_formula_hill = atoms.get_chemical_formula(mode='hill')
        except Exception as e:
            logger.warn('could not normalize material', exc_info=e)

        # Properties
        data = measurement.section_data
        if data:
            spectrum = data.section_spectrum
            if spectrum:
                if results.properties is None:
                    results.m_create(Properties)
                if results.properties.spectroscopy is None:
                    results.properties.m_create(SpectroscopyProperties)
                if results.method.method_name == 'EELS':
                    results.properties.spectroscopy.eels = spectrum
                else:
                    results.properties.spectroscopy.other_spectrum = spectrum

    def normalize_run(self, logger=None) -> None:
        # Fetch different information resources from which data is gathered
        repr_sys = None
        for section in self.section_run.system:
            if section.is_representative:
                repr_sys = section
                break
        try:
            encyclopedia = self.entry_archive.metadata.encyclopedia
        except Exception:
            encyclopedia = None
        try:
            optimade = self.entry_archive.metadata.optimade
        except Exception:
            optimade = None

        symmetry = None
        if repr_sys and repr_sys.symmetry:
            symmetry = repr_sys.symmetry[0]

        # Create the section and populate the subsections
        results = self.entry_archive.results
        results.material = self.material(repr_sys, symmetry, encyclopedia, optimade)
        results.method = self.method(encyclopedia)
        results.properties = self.properties(repr_sys, symmetry, encyclopedia)
        self.geometry_optimization(results.method, results.properties)

    def material(
            self,
            repr_sys: System,
            symmetry: SystemSymmetry,
            encyclopedia: EncyclopediaMetadata,
            optimade: OptimadeEntry) -> Material:
        """Returns a populated Material subsection."""
        material = Material()

        if repr_sys:
            material.structural_type = repr_sys.type
            names, counts = atomutils.get_hill_decomposition(repr_sys.atoms.labels, reduced=True)
            material.chemical_formula_reduced_fragments = [
                "{}{}".format(n, int(c) if c != 1 else "") for n, c in zip(names, counts)
            ]
        if encyclopedia:
            material.material_id = encyclopedia.material.material_id
            classes = encyclopedia.material.material_classification
            if classes:
                classifications = json.loads(classes)
                material.functional_type = classifications.get("material_class_springer")
                material.compound_type = classifications.get("compound_class_springer")
            material.material_name = encyclopedia.material.material_name
        if optimade:
            material.elements = optimade.elements
            material.chemical_formula_descriptive = optimade.chemical_formula_descriptive
            material.chemical_formula_reduced = optimade.chemical_formula_reduced
            material.chemical_formula_hill = optimade.chemical_formula_hill
            material.chemical_formula_anonymous = optimade.chemical_formula_anonymous

        if symmetry:
            symm = self.symmetry(repr_sys, symmetry, encyclopedia)
            if symm:
                material.symmetry = symm

        return material

    def symmetry(
            self,
            repr_sys: System,
            symmetry: SystemSymmetry,
            encyclopedia: EncyclopediaMetadata) -> Symmetry:
        """Returns a populated Symmetry subsection."""
        result = Symmetry()
        filled = False

        if symmetry:
            result.hall_number = symmetry.hall_number
            result.hall_symbol = symmetry.hall_symbol
            result.bravais_lattice = symmetry.bravais_lattice
            result.crystal_system = symmetry.crystal_system
            result.space_group_number = symmetry.space_group_number
            result.space_group_symbol = symmetry.international_short_symbol
            result.point_group = symmetry.point_group
            filled = True

        if encyclopedia and encyclopedia.material.bulk:
            result.strukturbericht_designation = encyclopedia.material.bulk.strukturbericht_designation
            result.structure_name = encyclopedia.material.bulk.structure_type
            result.prototype_formula = encyclopedia.material.bulk.structure_prototype
            filled = True

        proto = repr_sys.prototype if repr_sys else None
        proto = proto[0] if proto else None
        if proto:
            result.prototype_aflow_id = proto.aflow_id
            filled = True

        if filled:
            return result
        return None

    def wyckoff_sets(self, struct: StructureConventional, wyckoff_sets: Dict) -> None:
        """Populates the Wyckoff sets in the given structure.
        """
        for group in wyckoff_sets:
            wset = struct.m_create(WyckoffSet)
            if group.x is not None or group.y is not None or group.z is not None:
                if group.x is not None:
                    wset.x = float(group.x)
                if group.y is not None:
                    wset.y = float(group.y)
                if group.z is not None:
                    wset.z = float(group.z)
            wset.indices = group.indices
            wset.element = group.element
            wset.wyckoff_letter = group.wyckoff_letter

    def species(self, labels: List[str], struct: Structure) -> None:
        """Given a list of species labels, creates the corresponding Species
        sections in the given structure.
        """
        if labels is None:
            return
        unique_labels = sorted(list(set(labels)))
        for label in unique_labels:
            i_species = struct.m_create(Species)
            i_species.name = label
            try:
                symbol = label_to_chemical_symbol(label)
            except ValueError:
                self.logger.info("could not identify chemical symbol from the label: {}".format(label))
            else:
                i_species.chemical_symbols = [symbol]
            i_species.concentration = [1.0]

    def basis_set_type(self, repr_method) -> Union[str, None]:
        try:
            name = repr_method.basis_set[0].type
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

    def basis_set_name(self, repr_method) -> Union[str, None]:
        try:
            name = repr_method.basis_set[0].name
        except Exception:
            name = None
        return name

    def core_electron_treatment(self) -> str:
        treatment = config.services.unavailable_value
        code_name = self.section_run.program.name
        if code_name is not None:
            core_electron_treatments = {
                'VASP': 'pseudopotential',
                'FHI-aims': 'full all electron',
                'exciting': 'full all electron',
                'quantum espresso': 'pseudopotential'
            }
            treatment = core_electron_treatments.get(code_name, config.services.unavailable_value)
        return treatment

    def xc_functional_names(self, repr_method) -> Union[List[str], None]:
        if repr_method:
            functionals = []
            try:
                for functional_type in ['exchange', 'correlation', 'hybrid', 'contributions']:
                    functionals.extend([f.name for f in repr_method.dft.xc_functional[functional_type]])
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

    def method(
            self,
            encyclopedia: EncyclopediaMetadata) -> Method:
        """Returns a populated Method subsection."""
        method = Method()
        simulation = Simulation()
        repr_method = None
        method_name = config.services.unavailable_value
        methods = self.section_run.method
        n_methods = len(methods)

        def get_method_name(section_method):
            method_name = config.services.unavailable_value
            if section_method.electronic and section_method.electronic.method:
                method_name = section_method.electronic.method
            else:
                if section_method.gw is not None:
                    method_name = "GW"
            return method_name

        # If only one method is speficied use it directly
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

        if method_name == "GW":
            method.method_name = "GW"
            gw = GW()
            gw.type = repr_method.gw.type
            gw.starting_point = repr_method.gw.starting_point.split()
            simulation.gw = gw
        elif method_name in {"DFT", "DFT+U"}:
            method.method_name = "DFT"
            dft = DFT()
            dft.basis_set_type = self.basis_set_type(repr_method)
            dft.basis_set_name = self.basis_set_name(repr_method)
            dft.core_electron_treatment = self.core_electron_treatment()
            if repr_method.electronic is not None:
                if repr_method.electronic.smearing is not None:
                    dft.smearing_kind = repr_method.electronic.smearing.kind
                    dft.smearing_width = repr_method.electronic.smearing.width
                if repr_method.electronic.n_spin_channels:
                    dft.spin_polarized = repr_method.electronic.n_spin_channels > 1
                dft.van_der_Waals_method = repr_method.electronic.van_der_waals_method
                dft.relativity_method = repr_method.electronic.relativity_method
            dft.xc_functional_names = self.xc_functional_names(repr_method)
            dft.xc_functional_type = self.xc_functional_type(dft.xc_functional_names)
            if repr_method.scf is not None:
                dft.scf_threshold_energy_change = repr_method.scf.threshold_energy_change
            simulation.dft = dft

        if encyclopedia and encyclopedia.method:
            method.method_id = encyclopedia.method.method_id
        simulation.program_name = self.section_run.program.name
        simulation.program_version = self.section_run.program.version
        method.simulation = simulation

        return method

    def band_structure_electronic(self) -> Union[BandStructureElectronic, None]:
        """Returns a new section containing an electronic band structure. In
        the case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of kpoints.
          - There is a non-empty array of energies.
        """
        path = ["run", "calculation", "band_structure_electronic"]
        for bs in self.traverse_reversed(path):
            if not bs.segment:
                continue
            valid = True
            for segment in bs.segment:
                energies = segment.energies
                k_points = segment.kpoints
                if not valid_array(energies) or not valid_array(k_points):
                    valid = False
                    break
            if valid:
                # Fill band structure data to the newer, improved data layout
                bs_new = BandStructureElectronic()
                bs_new.reciprocal_cell = bs
                bs_new.segment = bs.segment
                bs_new.spin_polarized = bs_new.segment[0].energies.shape[0] > 1
                bs_new.energy_fermi = bs.energy_fermi
                for info in bs.band_gap:
                    info_new = bs_new.m_create(BandGap)
                    info_new.index = info.index
                    info_new.value = info.value
                    info_new.type = info.type
                    info_new.energy_highest_occupied = info.energy_highest_occupied
                    info_new.energy_lowest_unoccupied = info.energy_lowest_unoccupied
                return bs_new

        return None

    def dos_electronic(self) -> Union[DOSElectronic, None]:
        """Returns a reference to the section containing an electronic dos. In
        the case of multiple valid DOSes, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of dos_values_normalized.
          - There is a non-empty array of dos_energies.
        """
        path = ["run", "calculation", "dos_electronic"]
        for dos in self.traverse_reversed(path):
            energies = dos.energies
            values = np.array([d.value.magnitude for d in dos.total])
            if valid_array(energies) and valid_array(values):
                dos_new = DOSElectronic()
                dos_new.energies = dos
                dos_new.total = dos.total
                n_channels = values.shape[0]
                dos_new.spin_polarized = n_channels > 1
                dos_new.energy_fermi = dos.energy_fermi
                for info in dos.band_gap:
                    info_new = dos_new.m_create(BandGap)
                    info_new.index = info.index
                    info_new.energy_highest_occupied = info.energy_highest_occupied
                    info_new.energy_lowest_unoccupied = info.energy_lowest_unoccupied
                return dos_new

        return None

    def band_structure_phonon(self) -> Union[BandStructurePhonon, None]:
        """Returns a new section containing a phonon band structure. In
        the case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of kpoints.
          - There is a non-empty array of energies.
        """
        path = ["run", "calculation", "band_structure_phonon"]
        for bs in self.traverse_reversed(path):
            if not bs.segment:
                continue
            valid = True
            for segment in bs.segment:
                energies = segment.energies
                k_points = segment.kpoints
                if not valid_array(energies) or not valid_array(k_points):
                    valid = False
                    break
            if valid:
                # Fill band structure data to the newer, improved data layout
                bs_new = BandStructurePhonon()
                bs_new.segment = bs.segment
                return bs_new

        return None

    def dos_phonon(self) -> Union[DOSPhonon, None]:
        """Returns a section containing phonon dos data. In the case of
        multiple valid data sources, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of values.
          - There is a non-empty array of energies.
        """
        path = ["run", "calculation", "dos_phonon"]
        for dos in self.traverse_reversed(path):
            energies = dos.energies
            values = np.array([d.value.magnitude for d in dos.total])
            if valid_array(energies) and valid_array(values):
                dos_new = DOSPhonon()
                dos_new.energies = dos
                dos_new.total = dos.total
                return dos_new

        return None

    def energy_free_helmholtz(self) -> Union[EnergyFreeHelmholtz, None]:
        """Returns a section Helmholtz free energy data. In the case of
        multiple valid data sources, only the latest one is reported.

       Helmholtz free energy is reported only under the following conditions:
          - There is a non-empty array of temperatures.
          - There is a non-empty array of energies.
        """
        path = ["workflow", "thermodynamics"]
        for thermo_prop in self.traverse_reversed(path):
            temperatures = thermo_prop.temperature
            energies = thermo_prop.vibrational_free_energy_at_constant_volume
            if valid_array(temperatures) and valid_array(energies):
                energy_free = EnergyFreeHelmholtz()
                energy_free.energies = thermo_prop
                energy_free.temperatures = thermo_prop
                return energy_free

        return None

    def heat_capacity_constant_volume(self) -> Union[HeatCapacityConstantVolume, None]:
        """Returns a section containing heat capacity data. In the case of
        multiple valid data sources, only the latest one is reported.

       Heat capacity is reported only under the following conditions:
          - There is a non-empty array of temperatures.
          - There is a non-empty array of energies.
        """
        path = ["workflow", "thermodynamics"]
        for thermo_prop in self.traverse_reversed(path):
            temperatures = thermo_prop.temperature
            heat_capacities = thermo_prop.heat_capacity_c_v
            if valid_array(temperatures) and valid_array(heat_capacities):
                heat_cap = HeatCapacityConstantVolume()
                heat_cap.heat_capacities = thermo_prop
                heat_cap.temperatures = thermo_prop
                return heat_cap

        return None

    def geometry_optimization(self, method: Method, properties: Properties) -> None:
        """Populates both geometry optimization methodology and calculated
        properties based on the first found geometry optimization workflow.
        """
        path = ["workflow"]
        for workflow in self.traverse_reversed(path):
            # Check validity
            if workflow.type == "geometry_optimization" and workflow.calculations_ref:

                # Method
                geo_opt_wf = workflow.geometry_optimization
                if geo_opt_wf is not None:
                    geo_opt_meth = GeometryOptimizationMethod()
                    geo_opt_meth.type = geo_opt_wf.type
                    geo_opt_meth.convergence_tolerance_energy_difference = geo_opt_wf.convergence_tolerance_energy_difference
                    geo_opt_meth.convergence_tolerance_force_maximum = geo_opt_wf.convergence_tolerance_force_maximum
                    method.simulation.geometry_optimization = geo_opt_meth

                # Properties
                geo_opt_prop = GeometryOptimizationProperties()
                geo_opt_prop.trajectory = workflow.calculations_ref
                system_ref = workflow.calculation_result_ref.system_ref
                structure_optimized = self.structure_optimized(system_ref)
                if structure_optimized:
                    geo_opt_prop.structure_optimized = structure_optimized
                if geo_opt_wf is not None:
                    if geo_opt_wf.energies is not None:
                        geo_opt_prop.energies = geo_opt_wf
                    geo_opt_prop.final_energy_difference = geo_opt_wf.final_energy_difference
                    geo_opt_prop.final_force_maximum = geo_opt_wf.final_force_maximum
                properties.geometry_optimization = geo_opt_prop

                return

    def properties(
            self,
            repr_sys: System,
            symmetry: SystemSymmetry,
            encyclopedia: EncyclopediaMetadata) -> Properties:
        """Returns a populated Properties subsection."""
        properties = Properties()

        # Structures
        struct_orig = self.structure_original(repr_sys)
        struct_prim = self.structure_primitive(symmetry)
        struct_conv = self.structure_conventional(symmetry)
        if struct_orig or struct_prim or struct_conv:
            structures = Structures()
            if struct_conv:
                structures.structure_conventional = struct_conv
            if struct_prim:
                structures.structure_primitive = struct_prim
            if struct_orig:
                structures.structure_original = struct_orig
            properties.structures = structures

        # Electronic
        bs_electronic = self.band_structure_electronic()
        dos_electronic = self.dos_electronic()
        if bs_electronic or dos_electronic:
            electronic = ElectronicProperties()
            if bs_electronic:
                electronic.band_structure_electronic = bs_electronic
            if dos_electronic:
                electronic.dos_electronic = dos_electronic
            properties.electronic = electronic

        # Vibrational
        bs_phonon = self.band_structure_phonon()
        dos_phonon = self.dos_phonon()
        energy_free = self.energy_free_helmholtz()
        heat_cap = self.heat_capacity_constant_volume()
        if bs_phonon or dos_phonon or energy_free or heat_cap:
            vibrational = VibrationalProperties()
            if dos_phonon:
                vibrational.dos_phonon = dos_phonon
            if bs_phonon:
                vibrational.band_structure_phonon = bs_phonon
            if energy_free:
                vibrational.energy_free_helmholtz = energy_free
            if heat_cap:
                vibrational.heat_capacity_constant_volume = heat_cap
            properties.vibrational = vibrational

        # Mechanical
        energy_volume_curves = self.energy_volume_curves()
        bulk_modulus = self.bulk_modulus()
        shear_modulus = self.shear_modulus()
        if energy_volume_curves or bulk_modulus or shear_modulus:
            mechanical = MechanicalProperties()
            for ev in energy_volume_curves:
                mechanical.m_add_sub_section(MechanicalProperties.energy_volume_curve, ev)
            for bm in bulk_modulus:
                mechanical.m_add_sub_section(MechanicalProperties.bulk_modulus, bm)
            for sm in shear_modulus:
                mechanical.m_add_sub_section(MechanicalProperties.shear_modulus, sm)
            properties.mechanical = mechanical

        try:
            n_calc = len(self.section_run.calculation)
        except Exception:
            n_calc = 0
        properties.n_calculations = n_calc

        return properties

    def structure_original(self, repr_sys: System) -> StructureOriginal:
        """Returns a populated Structure subsection for the original
        structure.
        """
        if repr_sys and repr_sys.atoms:
            struct = StructureOriginal()
            struct.cartesian_site_positions = repr_sys.atoms.positions
            struct.species_at_sites = repr_sys.atoms.labels
            self.species(struct.species_at_sites, struct)
            lattice_vectors = repr_sys.atoms.lattice_vectors
            if lattice_vectors is not None and atomutils.is_valid_basis(lattice_vectors.magnitude):
                struct.dimension_types = np.array(repr_sys.atoms.periodic).astype(int)
                struct.lattice_vectors = lattice_vectors
                struct.cell_volume = atomutils.get_volume(lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(lattice_vectors)
            return struct

        return None

    def structure_primitive(self, symmetry: SystemSymmetry) -> StructurePrimitive:
        """Returns a populated Structure subsection for the primitive
        structure.
        """
        if symmetry:
            struct = StructurePrimitive()
            prim_sys = symmetry.system_primitive[0]
            struct.species_at_sites = atomutils.chemical_symbols(prim_sys.atomic_numbers)
            self.species(struct.species_at_sites, struct)
            lattice_vectors = prim_sys.lattice_vectors
            if lattice_vectors is not None and atomutils.is_valid_basis(lattice_vectors.magnitude):
                struct.cartesian_site_positions = atomutils.to_cartesian(prim_sys.positions.magnitude, lattice_vectors.magnitude)
                struct.dimension_types = [1, 1, 1]
                struct.lattice_vectors = lattice_vectors
                struct.cell_volume = atomutils.get_volume(lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(lattice_vectors)
            return struct

        return None

    def structure_conventional(self, symmetry: SystemSymmetry) -> StructureConventional:
        """Returns a populated Structure subsection for the conventional
        structure.
        """
        if symmetry:
            struct = StructureConventional()
            conv_sys = symmetry.system_std[0]
            struct.species_at_sites = atomutils.chemical_symbols(conv_sys.atomic_numbers)
            self.species(struct.species_at_sites, struct)
            lattice_vectors = conv_sys.lattice_vectors
            if lattice_vectors is not None and atomutils.is_valid_basis(lattice_vectors.magnitude):
                struct.cartesian_site_positions = atomutils.to_cartesian(conv_sys.positions.magnitude, lattice_vectors.magnitude)
                struct.dimension_types = [1, 1, 1]
                struct.lattice_vectors = lattice_vectors
                struct.cell_volume = atomutils.get_volume(lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(lattice_vectors)
                analyzer = symmetry.m_cache["symmetry_analyzer"]
                sets = analyzer.get_wyckoff_sets_conventional(return_parameters=True)
                self.wyckoff_sets(struct, sets)
            return struct

        return None

    def structure_optimized(self, system: System) -> Structure:
        """Returns a populated Structure subsection for the optimized
        structure.
        """
        if system:
            struct = StructureOptimized()
            struct.cartesian_site_positions = system.atoms.positions
            struct.species_at_sites = system.atoms.labels
            self.species(struct.species_at_sites, struct)
            lattice_vectors = system.atoms.lattice_vectors
            if lattice_vectors is not None and atomutils.is_valid_basis(lattice_vectors.magnitude):
                struct.dimension_types = np.array(system.atoms.periodic).astype(int)
                struct.lattice_vectors = system.atoms.lattice_vectors
                struct.cell_volume = atomutils.get_volume(system.atoms.lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(system.atoms.lattice_vectors)
            return struct

        return None

    def lattice_parameters(self, lattice_vectors) -> LatticeParameters:
        """Converts the given cell into LatticeParameters"""
        param_values = atomutils.cell_to_cellpar(lattice_vectors.magnitude)
        params = LatticeParameters()
        params.a = float(param_values[0])
        params.b = float(param_values[1])
        params.c = float(param_values[2])
        params.alpha = float(param_values[3])
        params.beta = float(param_values[4])
        params.gamma = float(param_values[5])
        return params

    def energy_volume_curves(self) -> List[EnergyVolumeCurve]:
        """Returns a list containing the found EnergyVolumeCurves.
        """
        workflows = self.entry_archive.workflow
        ev_curves = []
        for workflow in workflows:
            # Equation of state must be present
            equation_of_state = workflow.equation_of_state
            if not equation_of_state:
                continue

            # Volumes must be present
            volumes = equation_of_state.volumes
            if not valid_array(volumes):
                self.logger.warning("missing eos volumes")
                continue

            # Raw EV curve
            energies_raw = equation_of_state.energies
            if valid_array(energies_raw):
                ev_curves.append(EnergyVolumeCurve(
                    type="raw",
                    volumes=equation_of_state,
                    energies_raw=equation_of_state,
                ))
            else:
                self.logger.warning("missing eos energies")

            # Fitted EV curves
            fits = equation_of_state.eos_fit
            if not fits:
                continue
            for fit in fits:
                energies_fitted = fit.fitted_energies
                function_name = fit.function_name
                if valid_array(energies_fitted):
                    ev_curves.append(EnergyVolumeCurve(
                        type=function_name,
                        volumes=equation_of_state,
                        energies_fit=fit,
                    ))

        return ev_curves

    def bulk_modulus(self) -> List[BulkModulus]:
        """Returns a list containing the found BulkModulus.
        """
        workflows = self.entry_archive.workflow
        bulk_modulus = []
        for workflow in workflows:
            # From elastic workflow
            elastic = workflow.elastic
            if elastic:
                bulk_modulus_vrh = elastic.bulk_modulus_hill
                if bulk_modulus_vrh:
                    bulk_modulus.append(BulkModulus(
                        type="voigt_reuss_hill_average",
                        value=bulk_modulus_vrh,
                    ))
                bulk_modulus_voigt = elastic.bulk_modulus_voigt
                if bulk_modulus_voigt:
                    bulk_modulus.append(BulkModulus(
                        type="voigt_average",
                        value=bulk_modulus_voigt,
                    ))
                bulk_modulus_reuss = elastic.bulk_modulus_reuss
                if bulk_modulus_reuss:
                    bulk_modulus.append(BulkModulus(
                        type="reuss_average",
                        value=bulk_modulus_reuss,
                    ))

            # From energy-volume curve fit
            equation_of_state = workflow.equation_of_state
            if equation_of_state:
                fits = equation_of_state.eos_fit
                if not fits:
                    continue
                for fit in fits:
                    modulus = fit.bulk_modulus
                    function_name = fit.function_name
                    if modulus is not None and function_name:
                        bulk_modulus.append(BulkModulus(
                            type=function_name,
                            value=modulus,
                        ))
                    else:
                        self.logger.warning("missing eos fitted energies and/or function name")

        return bulk_modulus

    def shear_modulus(self) -> List[ShearModulus]:
        """Returns a list containing the found ShearModulus.
        """
        workflows = self.entry_archive.workflow
        shear_modulus = []
        for workflow in workflows:
            # From elastic workflow
            elastic = workflow.elastic
            if elastic:
                shear_modulus_vrh = elastic.shear_modulus_hill
                if shear_modulus_vrh:
                    shear_modulus.append(ShearModulus(
                        type="voigt_reuss_hill_average",
                        value=shear_modulus_vrh,
                    ))
                shear_modulus_voigt = elastic.shear_modulus_voigt
                if shear_modulus_voigt:
                    shear_modulus.append(ShearModulus(
                        type="voigt_average",
                        value=shear_modulus_voigt,
                    ))
                shear_modulus_reuss = elastic.shear_modulus_reuss
                if shear_modulus_reuss:
                    shear_modulus.append(ShearModulus(
                        type="reuss_average",
                        value=shear_modulus_reuss,
                    ))

        return shear_modulus

    def traverse_reversed(self, path: List[str]) -> Any:
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
        for t in traverse(self.entry_archive, path, 0):
            if t is not None:
                yield t
