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

import re
from nomad.datamodel.metainfo.workflow import Workflow
import numpy as np
from typing import List, Union, Any, Set, Optional
import ase.data
from matid import SymmetryAnalyzer
import matid.geometry

from nomad import config
from nomad import atomutils
from nomad.normalizing.common import structure_from_ase_atoms, structure_from_nomad_atoms, structures_2d
from nomad.normalizing.normalizer import Normalizer
from nomad.normalizing.method import MethodNormalizer
from nomad.normalizing.material import MaterialNormalizer
from nomad.datamodel.optimade import Species
from nomad.datamodel.metainfo.simulation.system import System, Symmetry as SystemSymmetry
from nomad.datamodel.results import (
    BandGapElectronic,
    RadialDistributionFunction,
    MeanSquaredDisplacement,
    Results,
    Material,
    Method,
    GeometryOptimization,
    Trajectory,
    MolecularDynamics,
    Methodology,
    TemperatureDynamic,
    VolumeDynamic,
    PressureDynamic,
    EnergyDynamic,
    Properties,
    StructuralProperties,
    DynamicalProperties,
    Structures,
    Structure,
    EnergyVolumeCurve,
    BulkModulus,
    ShearModulus,
    MechanicalProperties,
    ElectronicProperties,
    VibrationalProperties,
    ThermodynamicProperties,
    BandStructureElectronic,
    BandStructurePhonon,
    DOSElectronic,
    DOSPhonon,
    EnergyFreeHelmholtz,
    HeatCapacityConstantVolume,
)

re_label = re.compile("^([a-zA-Z][a-zA-Z]?)[^a-zA-Z]*")
elements = set(ase.data.chemical_symbols)


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

        for measurement in self.entry_archive.measurement:
            self.normalize_measurement(measurement)

        # Add the list of available_properties: it is a selected subset of the
        # stored properties.
        available_property_names = {
            "results.properties.electronic.band_structure_electronic.band_gap": "electronic.band_structure_electronic.band_gap",
            "results.properties.electronic.band_structure_electronic": "band_structure_electronic",
            "results.properties.electronic.dos_electronic": "dos_electronic",
            "results.properties.vibrational.dos_phonon": "dos_phonon",
            "results.properties.vibrational.band_structure_phonon": "band_structure_phonon",
            "results.properties.vibrational.energy_free_helmholtz": "energy_free_helmholtz",
            "results.properties.vibrational.heat_capacity_constant_volume": "heat_capacity_constant_volume",
            "results.properties.thermodynamic.trajectory": "trajectory",
            "results.properties.structural.radial_distribution_function": "radial_distribution_function",
            "results.properties.dynamical.mean_squared_displacement": "mean_squared_displacement",
            "results.properties.geometry_optimization": "geometry_optimization",
            "results.properties.mechanical.bulk_modulus": "bulk_modulus",
            "results.properties.mechanical.shear_modulus": "shear_modulus",
            "results.properties.mechanical.energy_volume_curve": "energy_volume_curve",
            "results.properties.spectroscopy.eels": "eels",
        }
        available_properties: List[str] = []
        for path, shortcut in available_property_names.items():
            for _ in self.traverse_reversed(path.split('.')):
                available_properties.append(shortcut)
                break
        results.properties.available_properties = sorted(available_properties)

    def normalize_sample(self, sample) -> None:
        results = self.entry_archive.results

        if results.material is None:
            results.material = Material()
        material = results.material
        if sample.elements and len(sample.elements) > 0:
            material.elements = sample.elements
        else:
            # Try to guess elements from sample formula or name
            if sample.chemical_formula:
                try:
                    material.elements = list(set(ase.Atoms(sample.chemical_formula).get_chemical_symbols()))
                except Exception:
                    if sample.name:
                        try:
                            material.elements = list(set(ase.Atoms(sample.name).get_chemical_symbols()))
                        except Exception:
                            pass
        if sample.chemical_formula:
            material.chemical_formula_descriptive = sample.chemical_formula

        try:
            material.elements = material.elements if material.elements else []
            atoms = None
            if material.chemical_formula_descriptive:
                try:
                    atoms = ase.Atoms(material.chemical_formula_descriptive)
                except Exception as e:
                    self.logger.warn('could not normalize formula, using elements next', exc_info=e)

            if atoms is None:
                atoms = ase.Atoms(''.join(material.elements))

            if material.elements is None or len(material.elements) == 0:
                material.elements = atoms.get_chemical_symbols()

            results.material.chemical_formula_hill = atoms.get_chemical_formula(mode='hill')
            results.material.chemical_formula_reduced = atoms.get_chemical_formula(mode='reduce')
            results.material.chemical_formula_descriptive = results.material.chemical_formula_hill
        except Exception as e:
            self.logger.warn('could not normalize material', exc_info=e)

    def normalize_measurement(self, measurement) -> None:
        results = self.entry_archive.results

        # Method
        if results.method is None:
            results.method = Method(
                method_name=measurement.method_abbreviation)

        # Sample
        if results.material is None:
            results.material = Material(elements=[])
        if len(measurement.sample) > 0:
            self.normalize_sample(measurement.sample[0])

    def normalize_run(self, logger=None) -> None:
        # Fetch different information resources from which data is gathered
        repr_system = None
        for section in self.section_run.system:
            if section.is_representative:
                repr_system = section
                break
        try:
            optimade = self.entry_archive.metadata.optimade
        except Exception:
            optimade = None

        repr_symmetry = None
        if repr_system and repr_system.symmetry:
            repr_symmetry = repr_system.symmetry[0]

        # Create the section and populate the subsections
        results = self.entry_archive.results
        properties, conv_atoms, wyckoff_sets, spg_number = self.properties(repr_system, repr_symmetry)
        results.properties = properties
        results.material = MaterialNormalizer(
            self.entry_archive,
            repr_system,
            repr_symmetry,
            spg_number,
            conv_atoms,
            wyckoff_sets,
            properties,
            optimade,
            logger
        ).material()
        results.method = MethodNormalizer(self.entry_archive, repr_system, results.material, logger).method()

    def species(self, labels: List[str], atomic_numbers: List[int], struct: Structure) -> None:
        """Given a list of species labels, creates the corresponding Species
        sections in the given structure.
        """
        if labels is None or atomic_numbers is None:
            return
        species: Set[str] = set()
        for label, atomic_number in zip(labels, atomic_numbers):
            if label not in species:
                species.add(label)
                i_species = struct.m_create(Species)
                i_species.name = label
                try:
                    symbol = atomutils.chemical_symbols([atomic_number])[0]
                except ValueError:
                    self.logger.info("could not identify chemical symbol for atomic number {}".format(atomic_number))
                else:
                    i_species.chemical_symbols = [symbol]
                i_species.concentration = [1.0]

    def band_structure_electronic(self) -> Union[List[BandStructureElectronic], None]:
        """Returns a new section containing an electronic band structure. In
        the case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of kpoints.
          - There is a non-empty array of energies.
        """
        def resolve_band_structure(path):
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
                        info_new = bs_new.m_create(BandGapElectronic)
                        info_new.index = info.index
                        info_new.value = info.value
                        info_new.type = info.type
                        info_new.energy_highest_occupied = info.energy_highest_occupied
                        info_new.energy_lowest_unoccupied = info.energy_lowest_unoccupied
                    return bs_new
            return None

        try:
            workflow_type = self.entry_archive.workflow[-1].type
        except Exception:
            workflow_type = None

        if workflow_type == 'GW':
            band_structures = []
            for method in ['gw', 'dft']:
                band_structure = resolve_band_structure(["workflow", "gw", f"band_structure_{method}"])
                if band_structure:
                    name = method.upper()
                    band_structure.label = name
                    for band_gap in band_structure.band_gap:
                        band_gap.label = name
                    band_structures.append(band_structure)
            return band_structures
        else:
            band_structure = resolve_band_structure(["run", "calculation", "band_structure_electronic"])
            if band_structure:
                return [band_structure]
        return None

    def dos_electronic(self) -> Union[List[DOSElectronic], None]:
        """Returns a reference to the section containing an electronic dos. In
        the case of multiple valid DOSes, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of dos_values_normalized.
          - There is a non-empty array of dos_energies.
        """

        def resolve_dos(path):
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
                        info_new = dos_new.m_create(BandGapElectronic)
                        info_new.index = info.index
                        info_new.energy_highest_occupied = info.energy_highest_occupied
                        info_new.energy_lowest_unoccupied = info.energy_lowest_unoccupied
                    return dos_new
            return None

        try:
            workflow_type = self.entry_archive.workflow[-1].type
        except Exception:
            workflow_type = None

        if workflow_type == 'GW':
            doss = []
            for method in ['gw', 'dft']:
                dos = resolve_dos(["workflow", "gw", f"dos_{method}"])
                if dos:
                    dos.label = method.upper()
                    doss.append(dos)
            return doss
        else:
            dos = resolve_dos(["run", "calculation", "dos_electronic"])
            if dos:
                return [dos]
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

    def geometry_optimization(self) -> Union[GeometryOptimization, None]:
        """Populates both geometry optimization methodology and calculated
        properties based on the first found geometry optimization workflow.
        """
        path = ["workflow"]
        for workflow in self.traverse_reversed(path):
            # Check validity
            if workflow.type == "geometry_optimization" and workflow.calculations_ref:

                geo_opt = GeometryOptimization()
                geo_opt_wf = workflow.geometry_optimization
                geo_opt.trajectory = workflow.calculations_ref
                system_ref = workflow.calculation_result_ref.system_ref
                structure_optimized = structure_from_nomad_atoms(system_ref)
                if structure_optimized:
                    geo_opt.structure_optimized = structure_optimized
                if geo_opt_wf is not None:
                    geo_opt.type = geo_opt_wf.type
                    geo_opt.convergence_tolerance_energy_difference = geo_opt_wf.convergence_tolerance_energy_difference
                    geo_opt.convergence_tolerance_force_maximum = geo_opt_wf.convergence_tolerance_force_maximum
                    if geo_opt_wf.energies is not None:
                        geo_opt.energies = geo_opt_wf
                    geo_opt.final_energy_difference = geo_opt_wf.final_energy_difference
                    geo_opt.final_force_maximum = geo_opt_wf.final_force_maximum
                    geo_opt.final_displacement_maximum = geo_opt_wf.final_displacement_maximum
                return geo_opt

        return None

    def get_md_methodology(self, workflow: Workflow) -> Optional[MolecularDynamics]:
        """Retrieves the MD methodology from the given workflow.
        """
        md_wf = workflow.molecular_dynamics
        md = None
        if md_wf is not None:
            try:
                md = MolecularDynamics()
                md.time_step = md_wf.integration_parameters.integration_timestep
                md.ensemble_type = md_wf.thermodynamic_ensemble
            except Exception:
                pass
        return md

    def trajectory(self) -> List[Trajectory]:
        """Returns a list of trajectories.
        """
        path = ["workflow"]
        trajs = []
        for workflow in self.traverse_reversed(path):
            # Check validity
            if workflow.type == "molecular_dynamics":
                traj = Trajectory()
                md = self.get_md_methodology(workflow)
                if md:
                    traj.methodology = Methodology(molecular_dynamics=md)

                # Loop through calculations, gather thermodynamics directly
                # from each step in the workflow.
                volume = []
                volume_time = []
                pressure = []
                pressure_time = []
                temperature = []
                temperature_time = []
                potential_energy = []
                potential_energy_time = []

                for calc in workflow.calculations_ref:
                    time = calc.time
                    if time is not None:
                        time = time.magnitude
                        if calc.volume is not None:
                            volume.append(calc.volume.magnitude)
                            volume_time.append(time)
                        if calc.pressure is not None:
                            pressure.append(calc.pressure.magnitude)
                            pressure_time.append(time)
                        if calc.temperature is not None:
                            temperature.append(calc.temperature.magnitude)
                            temperature_time.append(time)
                        if calc.energy:
                            if calc.energy.potential is not None:
                                potential_energy.append(calc.energy.potential.value.magnitude)
                                potential_energy_time.append(time)

                available_properties = []
                if volume:
                    traj.volume = VolumeDynamic(value=volume, time=volume_time)
                    available_properties.append('volume')
                if pressure:
                    traj.pressure = PressureDynamic(value=pressure, time=pressure_time)
                    available_properties.append('pressure')
                if temperature:
                    traj.temperature = TemperatureDynamic(value=temperature, time=temperature_time)
                    available_properties.append('temperature')
                if potential_energy:
                    traj.energy_potential = EnergyDynamic(value=potential_energy, time=potential_energy_time)
                    available_properties.append('energy_potential')
                if available_properties:
                    traj.available_properties = available_properties
                trajs.append(traj)
        return trajs

    def rdf(self) -> List[RadialDistributionFunction]:
        """Returns a list of radial distribution functions.
        """
        path = ["workflow", "molecular_dynamics", "results", "radial_distribution_functions"]
        rdfs = []
        for rdf_workflow in self.traverse_reversed(path):
            rdf_values = rdf_workflow.radial_distribution_function_values
            if rdf_values is not None:
                for rdf_value in rdf_values or []:
                    rdf = RadialDistributionFunction()
                    try:
                        rdf.bins = rdf_value.bins
                        rdf.n_bins = rdf_value.n_bins
                        rdf.value = rdf_value.value
                        rdf.label = rdf_value.label
                        rdf.frame_start = rdf_value.frame_start
                        rdf.frame_end = rdf_value.frame_end
                        rdf.type = rdf_workflow.type
                        md = self.get_md_methodology(rdf_workflow.m_parent.m_parent.m_parent)
                        if md:
                            rdf.methodology = Methodology(
                                molecular_dynamics=md
                            )
                    except Exception as e:
                        self.logger.error('error in resolving radial distribution data', exc_info=e)
                    else:
                        rdfs.append(rdf)

        return rdfs

    def msd(self) -> List[MeanSquaredDisplacement]:
        """Returns a list of mean squared displacements.
        """
        path = ["workflow", "molecular_dynamics", "results", "mean_squared_displacements"]
        msds = []
        for msd_workflow in self.traverse_reversed(path):
            msd_values = msd_workflow.mean_squared_displacement_values
            if msd_values is not None:
                for msd_value in msd_values or []:
                    msd = MeanSquaredDisplacement()
                    try:
                        msd.times = msd_value.times
                        msd.n_times = msd_value.n_times
                        msd.value = msd_value.value
                        msd.label = msd_value.label
                        msd.errors = msd_value.errors
                        msd.type = msd_workflow.type
                        msd.direction = msd_workflow.direction
                        msd.error_type = msd_workflow.error_type
                        diffusion_constant = msd_value.diffusion_constant
                        if diffusion_constant is not None:
                            msd.diffusion_constant_value = diffusion_constant.value
                            msd.diffusion_constant_error_type = diffusion_constant.error_type
                            msd.diffusion_constant_errors = diffusion_constant.errors

                        md = self.get_md_methodology(msd_workflow.m_parent.m_parent.m_parent)
                        if md:
                            msd.methodology = Methodology(
                                molecular_dynamics=md
                            )
                    except Exception as e:
                        self.logger.error('error in resolving mean squared displacement data', exc_info=e)
                    else:
                        msds.append(msd)

        return msds

    def properties(
            self,
            repr_system: System,
            repr_symmetry: SystemSymmetry) -> tuple:
        """Returns a populated Properties subsection."""
        properties = Properties()

        # Structures
        struct_orig = None
        struct_prim = None
        struct_conv = None
        conv_atoms = None
        wyckoff_sets = None
        spg_number = None
        if repr_system:
            original_atoms = repr_system.m_cache.get("representative_atoms")
            if original_atoms:
                prim_atoms = None
                structural_type = repr_system.type
                if structural_type == "bulk":
                    conv_atoms, prim_atoms, wyckoff_sets, spg_number = self.structures_bulk(repr_symmetry)
                elif structural_type == "2D":
                    conv_atoms, prim_atoms, wyckoff_sets, spg_number = structures_2d(original_atoms)
                elif structural_type == "1D":
                    conv_atoms, prim_atoms = self.structures_1d(original_atoms)

                struct_orig = structure_from_ase_atoms(original_atoms, logger=self.logger)
                struct_prim = structure_from_ase_atoms(prim_atoms, logger=self.logger)
                wyckoff_sets_serialized = wyckoff_sets if structural_type == "bulk" else None
                struct_conv = structure_from_ase_atoms(conv_atoms, wyckoff_sets_serialized, logger=self.logger)

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
        geometry_optimization = self.geometry_optimization()
        if energy_volume_curves or bulk_modulus or shear_modulus or geometry_optimization:
            mechanical = MechanicalProperties()
            for ev in energy_volume_curves:
                mechanical.m_add_sub_section(MechanicalProperties.energy_volume_curve, ev)
            for bm in bulk_modulus:
                mechanical.m_add_sub_section(MechanicalProperties.bulk_modulus, bm)
            for sm in shear_modulus:
                mechanical.m_add_sub_section(MechanicalProperties.shear_modulus, sm)
            properties.mechanical = mechanical

        # Geometry optimization
        properties.geometry_optimization = self.geometry_optimization()

        # Thermodynamic
        trajectory = self.trajectory()
        if trajectory:
            thermodynamic = ThermodynamicProperties()
            thermodynamic.trajectory = trajectory
            properties.thermodynamic = thermodynamic

        # Structural
        rdf = self.rdf()
        if rdf:
            structural = StructuralProperties()
            structural.radial_distribution_function = rdf
            properties.structural = structural

        # Dynamical
        msd = self.msd()
        if msd:
            dynamical = DynamicalProperties()
            dynamical.mean_squared_displacement = msd
            properties.dynamical = dynamical

        try:
            n_calc = len(self.section_run.calculation)
        except Exception:
            n_calc = 0
        properties.n_calculations = n_calc

        return properties, conv_atoms, wyckoff_sets, spg_number

    def structures_bulk(self, repr_symmetry):
        """The symmetry of bulk structures has already been analyzed. Here we
        use the cached results.
        """
        conv_atoms = None
        prim_atoms = None
        wyckoff_sets = None
        spg_number = None
        if repr_symmetry:
            symmetry_analyzer = repr_symmetry.m_cache.get("symmetry_analyzer")
            if symmetry_analyzer:
                spg_number = symmetry_analyzer.get_space_group_number()
                conv_atoms = symmetry_analyzer.get_conventional_system()
                prim_atoms = symmetry_analyzer.get_primitive_system()

                # For some reason MatID seems to drop the periodicity, reintroduce it here.
                conv_atoms.set_pbc(True)
                prim_atoms.set_pbc(True)
                try:
                    wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional(return_parameters=True)
                except Exception:
                    self.logger.error('Error resolving Wyckoff sets.')
                    wyckoff_sets = []

        return conv_atoms, prim_atoms, wyckoff_sets, spg_number

    def structures_1d(self, original_atoms):
        conv_atoms = None
        prim_atoms = None
        try:
            # First get a symmetry analyzer and the primitive system
            symm_system = original_atoms.copy()
            symm_system.set_pbc(True)
            symmetry_analyzer = SymmetryAnalyzer(
                symm_system,
                config.normalize.symmetry_tolerance,
                config.normalize.flat_dim_threshold
            )
            prim_atoms = symmetry_analyzer.get_primitive_system()
            prim_atoms.set_pbc(True)

            # Get dimension of system by also taking into account the covalent radii
            dimensions = matid.geometry.get_dimensions(prim_atoms, [True, True, True])
            basis_dimensions = np.linalg.norm(prim_atoms.get_cell(), axis=1)
            gaps = basis_dimensions - dimensions
            periodicity = gaps <= config.normalize.cluster_threshold

            # If one axis is not periodic, return. This only happens if the vacuum
            # gap is not aligned with a cell vector.
            if sum(periodicity) != 1:
                self.logger.error("could not detect the periodic dimensions in a 1D system")
                return conv_atoms, prim_atoms

            # Translate to center of mass
            conv_atoms = prim_atoms.copy()
            pbc_cm = matid.geometry.get_center_of_mass(prim_atoms)
            cell_center = 0.5 * np.sum(conv_atoms.get_cell(), axis=0)
            translation = cell_center - pbc_cm
            translation[periodicity] = 0
            conv_atoms.translate(translation)
            conv_atoms.wrap()
            conv_atoms.set_pbc(periodicity)

            # Reduce cell size to just fit the system in the non-periodic dimensions.
            conv_atoms = atomutils.get_minimized_structure(conv_atoms)

            # Swap the cell axes so that the periodic one is always the first
            # basis (=a)
            swap_dim = 0
            for i, periodic in enumerate(periodicity):
                if periodic:
                    periodic_dim = i
                    break
            if periodic_dim != swap_dim:
                atomutils.swap_basis(conv_atoms, periodic_dim, swap_dim)

            prim_atoms = conv_atoms
        except Exception as e:
            self.logger.error(
                'could not construct a conventional system for a 1D material',
                exc_info=e
            )
        return conv_atoms, prim_atoms

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
            if not root:
                return
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
