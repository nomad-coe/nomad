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
import numpy as np
from ase import Atoms
from nptyping import NDArray

from nomad.datamodel import ArchiveSection, EntryArchive
from nomad.metainfo import MSection, SubSection, Section, Quantity, MEnum, Reference, derived
from nomad.datamodel.metainfo.common import FastAccess
from nomad.datamodel.metainfo.workflow2 import Workflow, Link, Task
from nomad.datamodel.metainfo.simulation.system import System, AtomsGroup
from nomad.datamodel.metainfo.simulation.method import (
    Method, XCFunctional, BasisSet, GW as GWMethodology
)
from nomad.datamodel.metainfo.simulation.calculation import (
    Calculation, Dos, BandStructure, BandEnergies, Density, Potential, Spectra,
    RadiusOfGyration as RadiusOfGyrationCalculation,
    RadiusOfGyrationValues as RadiusOfGyrationValuesCalculation)
from nomad.atomutils import archive_to_universe
from nomad.atomutils import (
    calc_molecular_rdf,
    calc_molecular_mean_squared_displacements,
    calc_molecular_radius_of_gyration)


def resolve_difference(values):
    delta_values = None

    for n in range(-1, -len(values), -1):
        a = values[n]
        b = values[n - 1]
        if a is None or b is None:
            continue
        delta_values = abs(a - b)
        if delta_values != 0.0:
            break

    return delta_values


_input_structure_name = 'Input structure'
_output_structure_name = 'Output structure'
_input_method_name = 'Input method'
_output_calculation_name = 'Output calculation'
_workflow_method_name = 'Workflow parameters'
_workflow_results_name = 'Workflow results'


class SimulationWorkflowMethod(ArchiveSection):

    m_def = Section(validate=False)

    def normalize(self, archive, logger):
        pass


class SimulationWorkflowResults(ArchiveSection):

    m_def = Section(validate=False)

    calculation_result_ref = Quantity(
        type=Reference(Calculation.m_def),
        shape=[],
        description='''
        Reference to calculation result. In the case of serial workflows, this corresponds
        to the final step in the simulation. For the parallel case, it refers to the original system.
        ''',
        categories=[FastAccess])

    def normalize(self, archive, logger):
        pass


class SimulationWorkflow(Workflow):

    method = SubSection(sub_section=SimulationWorkflowMethod)

    results = SubSection(sub_section=SimulationWorkflowResults)

    def normalize(self, archive: EntryArchive, logger):
        self._calculations: List[Calculation] = []
        self._systems: List[System] = []
        self._methods: List[Method] = []
        try:
            self._calculations = archive.run[-1].calculation
            self._systems = archive.run[-1].system
            self._methods = archive.run[-1].method
        except Exception:
            logger.warning('System, method and calculation required for normalization.')
            pass

        if not self._calculations or not self._systems:
            return

        if not self.inputs:
            if self._systems:
                self.m_add_sub_section(
                    Workflow.inputs, Link(name=_input_structure_name, section=self._systems[0]))

            if self.method:
                self.m_add_sub_section(
                    Workflow.inputs, Link(name=_workflow_method_name, section=self.method))

        for link in self.inputs:
            if isinstance(link.section, System):
                self.input_structure = link.section
                break

        if not self.outputs:
            if self._calculations:
                self.m_add_sub_section(
                    Workflow.outputs, Link(name=_output_calculation_name, section=self._calculations[-1]))

            if self.results:
                self.m_add_sub_section(
                    Workflow.outputs, Link(name=_workflow_results_name, section=self.results))


class Decomposition(MSection):
    '''
    Section containing information about the system to which an unstable compound will
    decompose to.
    '''

    m_def = Section(validate=False)

    fraction = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Amount of the resulting system.
        ''')

    system_ref = Quantity(
        type=Reference(System.m_def),
        shape=[],
        description='''
        Reference to the resulting system.
        ''')

    formula = Quantity(
        type=str,
        shape=[],
        description='''
        Chemical formula of the resulting system.
        ''')


class Stability(MSection):
    '''
    Section containing information regarding the stability of the system.
    '''

    m_def = Section(validate=False)

    n_references = Quantity(
        type=int,
        shape=[],
        description='''
        Number of reference systems.
        ''')

    systems_ref = Quantity(
        type=Reference(System.m_def),
        shape=['n_references'],
        description='''
        References to the reference systems.
        ''')

    formation_energy = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Calculated value of the formation energy of the compound.
        ''')

    delta_formation_energy = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Energy with respect to the convex hull.
        ''')

    n_references = Quantity(
        type=int,
        shape=[],
        description='''
        Number of reference systems.
        ''')

    is_stable = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if a compound is stable.
        ''')

    decomposition = SubSection(sub_section=Decomposition.m_def, repeats=True)


class ThermodynamicsResults(SimulationWorkflowResults):

    n_values = Quantity(
        type=int,
        shape=[],
        description='''
        Number of thermodynamics property evaluations.
        ''')

    temperature = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='kelvin',
        description='''
        Specifies the temperatures at which properties such as the Helmholtz free energy
        are calculated.
        ''')

    pressure = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='pascal',
        description='''
        Array containing the values of the pressure (one third of the trace of the stress
        tensor) corresponding to each property evaluation.
        ''')

    helmholtz_free_energy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Helmholtz free energy per unit cell at constant volume.
        ''')

    heat_capacity_c_p = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule / kelvin',
        description='''
        Heat capacity per cell unit at constant pressure.
        ''')

    heat_capacity_c_v = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule / kelvin',
        description='''
        Heat capacity per cell unit at constant volume.
        ''')

    @derived(
        type=np.float64,
        shape=['n_values'],
        unit='joule / (kelvin * kilogram)',
        description='''
        Specific heat capacity at constant volume.
        ''',
        cached=True
    )
    def heat_capacity_c_v_specific(self) -> NDArray:
        """Returns the specific heat capacity by dividing the heat capacity per
        cell with the mass of the atoms in the cell.
        """
        import nomad.atomutils
        workflow = self.m_parent
        if not workflow._systems or not workflow._systems[0].atoms:
            return
        atomic_numbers = workflow._systems[0].atoms.species
        mass_per_unit_cell = nomad.atomutils.get_summed_atomic_mass(atomic_numbers)
        heat_capacity = self.heat_capacity_c_v
        specific_heat_capacity = heat_capacity / mass_per_unit_cell

        return specific_heat_capacity.magnitude

    vibrational_free_energy_at_constant_volume = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Holds the vibrational free energy per cell unit at constant volume.
        ''')

    @derived(
        type=np.float64,
        shape=['n_values'],
        unit='joule / kilogram',
        description='''
        Stores the specific vibrational free energy at constant volume.
        ''',
        cached=True
    )
    def vibrational_free_energy_at_constant_volume_specific(self) -> NDArray:
        import nomad.atomutils
        workflow = self.m_parent
        if not workflow._systems or not workflow._systems[0].atoms:
            return
        atomic_numbers = workflow._systems[0].atoms.species
        mass_per_unit_cell = nomad.atomutils.get_summed_atomic_mass(atomic_numbers)
        free_energy = self.vibrational_free_energy_at_constant_volume
        specific_free_energy = free_energy / mass_per_unit_cell

        return specific_free_energy.magnitude

    vibrational_free_energy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Calculated value of the vibrational free energy, F_vib.
        ''')

    vibrational_internal_energy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Calculated value of the vibrational internal energy, U_vib.
        ''')

    vibrational_entropy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule / kelvin',
        description='''
        Calculated value of the vibrational entropy, S.
        ''')

    gibbs_free_energy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Calculated value of the Gibbs free energy, G.
        ''')

    entropy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule / kelvin',
        description='''
        Calculated value of the entropy.
        ''')

    enthalpy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Calculated value of enthalpy.
        ''')

    internal_energy = Quantity(
        type=np.float64,
        shape=['n_values'],
        unit='joule',
        description='''
        Calculated value of the internal energy, U.
        ''')

    stability = SubSection(sub_section=Stability.m_def, repeats=False)

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        try:
            calculations = archive.run[-1].calculation
        except Exception:
            calculations = []

        def set_thermo_property(name):
            values = []
            quantity = None
            for calc in calculations:
                if hasattr(calc, name):
                    try:
                        quantity = calc[name]
                        values.append(quantity.magnitude if hasattr(quantity, 'magnitude') else quantity)
                        continue
                    except Exception:
                        pass
                # TODO section thermodynamics should be removed
                for thermo in calc.thermodynamics:
                    try:
                        quantity = thermo[name]
                        values.append(quantity.magnitude if hasattr(quantity, 'magnitude') else quantity)
                    except Exception:
                        pass
            if len(values) == 0:
                return

            unit = quantity.units if hasattr(quantity, 'units') else 1.0
            setattr(self, name, np.array(values) * unit)

        if self.temperature is None:
            set_thermo_property('temperature')

        if self.helmholtz_free_energy is None:
            set_thermo_property('helmholtz_free_energy')

        if self.vibrational_free_energy_at_constant_volume is None:
            set_thermo_property('vibrational_free_energy_at_constant_volume')

        if self.heat_capacity_c_v is None:
            set_thermo_property('heat_capacity_c_v')


class SinglePointResults(SimulationWorkflowResults):

    n_scf_steps = Quantity(
        type=int,
        shape=[],
        description='''
        Number of self-consistent steps in the calculation.
        ''')

    final_scf_energy_difference = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        The difference in the energy between the last two scf steps.
        ''')

    is_converged = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the convergence criteria were fullfilled.
        ''')

    n_data = Quantity(
        type=np.int32,
        shape=[],
        description='''
        ''')

    dos = Quantity(
        type=Reference(Dos),
        shape=['n_data'],
        description='''
        Reference to the electronic density of states data.
        ''')

    band_structure = Quantity(
        type=Reference(BandStructure),
        shape=['n_data'],
        description='''
        Reference to the electronic band structure data.
        ''')

    eigenvalues = Quantity(
        type=Reference(BandEnergies),
        shape=['n_data'],
        description='''
        Reference to the eigenvalues.
        ''')

    potential = Quantity(
        type=Reference(Potential),
        shape=['n_data'],
        description='''
        Reference to the potential data.
        ''')

    density_charge = Quantity(
        type=Reference(Density),
        shape=['n_data'],
        description='''
        Reference to the charge density data.
        ''')

    spectra = Quantity(
        type=Reference(Spectra),
        shape=['n_data'],
        description='''
        Reference to the spectral data.
        ''')


class SinglePointMethod(SimulationWorkflowMethod):

    method = Quantity(
        type=str,
        shape=[],
        description='''
        Calculation method used.
        ''')


class SinglePoint(SimulationWorkflow):

    method = SubSection(sub_section=SinglePointMethod)

    results = SubSection(sub_section=SinglePointResults)

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.tasks:
            task = Task()
            if self._systems:
                task.m_add_sub_section(
                    Task.inputs, Link(name=_input_structure_name, section=self._systems[0]))
            if self._methods:
                task.m_add_sub_section(
                    Task.inputs, Link(name=_input_method_name, section=self._methods[0]))
            if self._calculations:
                task.m_add_sub_section(
                    Task.inputs, Link(name=_output_calculation_name, section=self._calculations[0]))

            self.tasks = [task]

        if not self.method:
            self.method = SinglePointMethod()

        if not self.inputs:
            self.m_add_sub_section(
                SimulationWorkflow.inputs, Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = SinglePointResults()

        if not self.outputs:
            self.m_add_sub_section(
                SimulationWorkflow.outputs, Link(name=_workflow_results_name, section=self.results))

        if not self.method.method:
            try:
                # TODO keep extending for other SinglePoint
                for method_name in ['dft', 'gw', 'bse', 'dmft']:
                    if self._methods[-1].m_xpath(method_name):
                        self.method.method = method_name.upper()
                        break
            except Exception:
                pass

        if not self._calculations:
            return

        last_calc = self._calculations[-1]
        if not self.results.n_scf_steps:
            self.results.n_scf_steps = len(last_calc.scf_iteration)

        energies = [scf.energy.total.value for scf in last_calc.scf_iteration if scf.energy is not None and scf.energy.total is not None]
        delta_energy = resolve_difference(energies)
        if not self.results.final_scf_energy_difference and delta_energy is not None:
            self.results.final_scf_energy_difference = delta_energy

        if not self.results.is_converged and delta_energy is not None:
            try:
                threshold = self._methods[-1].scf.threshold_energy_change
                self.results.is_converged = bool(delta_energy <= threshold)
            except Exception:
                pass

        if not self.results.dos and last_calc.dos_electronic:
            self.results.dos = last_calc.dos_electronic

        if not self.results.band_structure and last_calc.band_structure_electronic:
            self.results.band_structure = last_calc.band_structure_electronic

        if not self.results.eigenvalues and last_calc.eigenvalues:
            self.results.eigenvalues = last_calc.eigenvalues

        if not self.results.density_charge and last_calc.density_charge:
            self.results.density_charge = last_calc.density_charge

        if not self.results.potential and last_calc.potential:
            self.results.potential = last_calc.potential

        if not self.results.spectra and last_calc.spectra:
            self.results.spectra = last_calc.spectra

        if not self.results.calculation_result_ref:
            self.results.calculation_result_ref = last_calc


class ParallelSimulation(SimulationWorkflow):
    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.tasks:
            for n, calculation in enumerate(self._calculations):
                inputs, outputs = [], [Link(name=_output_calculation_name, section=calculation)]
                if self._calculations[n].system_ref:
                    inputs.append(Link(name=_input_structure_name, section=self._calculations[n].system_ref))
                elif len(self._calculations) == len(self._systems):
                    inputs.append(Link(name=_input_structure_name, section=self._systems[n]))
                else:
                    continue
                if self._calculations[n].method_ref:
                    inputs.append(Link(name=_input_method_name, section=self._calculations[n].method_ref))
                elif len(self._calculations) == len(self._methods):
                    inputs.append(Link(name=_input_method_name, section=self._methods[n]))
                elif len(self._methods) == 1:
                    inputs.append(Link(name=_input_method_name, section=self._methods[0]))
                self.tasks.append(Task(name=f'Calculation {n}', inputs=inputs, outputs=outputs))


class SerialSimulation(SimulationWorkflow):
    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if not self.tasks:
            previous_structure = None
            for n, calculation in enumerate(self._calculations):
                inputs, outputs = [], [Link(name=_output_calculation_name, section=calculation)]
                if calculation.system_ref:
                    input_structure = self.input_structure if n == 0 else self._calculations[n - 1].system_ref
                    if not input_structure:
                        input_structure = previous_structure
                    if input_structure:
                        inputs.append(Link(name=_input_structure_name, section=input_structure))
                    previous_structure = calculation.system_ref
                    outputs.append(Link(name=_output_structure_name, section=calculation.system_ref))
                elif len(self._calculations) == len(self._systems):
                    inputs.append(Link(name=_input_structure_name, section=self.input_structure if n == 0 else self._systems[n - 1]))
                    outputs.append(Link(name=_output_structure_name, section=self._systems[n]))
                else:
                    continue
                if calculation.method_ref:
                    inputs.append(Link(name=_input_method_name, section=calculation.method_ref))
                elif len(self._calculations) == len(self._methods):
                    inputs.append(Link(name=_input_method_name, section=self._methods[n]))
                elif len(self._methods) == 1:
                    inputs.append(Link(name=_input_method_name, section=self._methods[0]))
                self.tasks.append(Task(name=f'Step {n}', inputs=inputs, outputs=outputs))


class GeometryOptimizationMethod(SimulationWorkflowMethod):

    type = Quantity(
        type=MEnum('static', 'atomic', 'cell_shape', 'cell_volume'),
        shape=[],
        description='''
        The type of geometry optimization, which denotes what is being optimized.

        Allowed values are:

        | Type                   | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `"static"`             | no optimization |

        | `"atomic"`             | the atomic coordinates alone are updated |

        | `"cell_volume"`         | `"atomic"` + cell lattice paramters are updated isotropically |

        | `"cell_shape"`        | `"cell_volume"` but without the isotropic constraint: all cell parameters are updated |

        ''')

    method = Quantity(
        type=str,
        shape=[],
        description='''
        The method used for geometry optimization. Some known possible values are:
        `"steepest_descent"`, `"conjugant_gradient"`, `"low_memory_broyden_fletcher_goldfarb_shanno"`.
        ''')

    convergence_tolerance_energy_difference = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        The input energy difference tolerance criterion.
        ''')

    convergence_tolerance_force_maximum = Quantity(
        type=np.float64,
        shape=[],
        unit='newton',
        description='''
        The input maximum net force tolerance criterion.
        ''')

    convergence_tolerance_stress_maximum = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        The input maximum stress tolerance criterion.
        ''')

    convergence_tolerance_displacement_maximum = Quantity(
        type=np.float64,
        shape=[],
        unit='meter',
        description='''
        The input maximum displacement tolerance criterion.
        ''')

    optimization_steps_maximum = Quantity(
        type=int,
        shape=[],
        description='''
        Maximum number of optimization steps.
        ''')


class GeometryOptimizationResults(SimulationWorkflowResults):

    optimization_steps = Quantity(
        type=int,
        shape=[],
        description='''
        Number of saved optimization steps.
        ''')

    energies = Quantity(
        type=np.float64,
        unit='joule',
        shape=['optimization_steps'],
        description='''
        List of energy_total values gathered from the single configuration
        calculations that are a part of the optimization trajectory.
        ''')

    steps = Quantity(
        type=np.int32,
        shape=['optimization_steps'],
        description='''
        The step index corresponding to each saved configuration.
        ''')

    final_energy_difference = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        The difference in the energy_total between the last two steps during
        optimization.
        ''')

    final_force_maximum = Quantity(
        type=np.float64,
        shape=[],
        unit='newton',
        description='''
        The maximum net force in the last optimization step.
        ''')

    final_displacement_maximum = Quantity(
        type=np.float64,
        shape=[],
        unit='meter',
        description='''
        The maximum displacement in the last optimization step with respect to previous.
        ''')

    is_converged_geometry = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the geometry convergence criteria were fulfilled.
        ''')


class GeometryOptimization(SerialSimulation):

    method = SubSection(sub_section=GeometryOptimizationMethod)

    results = SubSection(sub_section=GeometryOptimizationResults)

    def _get_geometry_optimization_type(self):
        if not self._systems:
            return

        def compare_cell(cell1, cell2):
            if (cell1 == cell2).all():
                return None
            else:
                cell1_normed = cell1 / np.linalg.norm(cell1)
                cell2_normed = cell2 / np.linalg.norm(cell2)
                if (cell1_normed == cell2_normed).all():
                    return 'cell_volume'
                else:
                    return 'cell_shape'

        if len(self._systems) < 2:
            return 'static'

        else:
            if self._systems[0].atoms is None or self._systems[-1].atoms is None:
                return 'static'

            cell_init = self._systems[0].atoms.lattice_vectors
            cell_final = self._systems[-1].atoms.lattice_vectors
            if cell_init is None or cell_final is None:
                return 'static'

            cell_relaxation = compare_cell(cell_init.magnitude, cell_final.magnitude)

            if cell_relaxation is not None:
                return cell_relaxation

            atom_pos_init = self._systems[0].atoms.positions
            atom_pos_final = self._systems[-1].atoms.positions
            if atom_pos_init is None or atom_pos_final is None:
                return 'static'

            if (atom_pos_init.magnitude == atom_pos_final.magnitude).all():
                return 'static'

            return 'atomic'

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = GeometryOptimizationMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = GeometryOptimizationResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))

        if not self.method.type and self._systems:
            self.method.type = self._get_geometry_optimization_type()

        if not self.results.optimization_steps:
            self.results.optimization_steps = len(self._calculations)

        energies = []
        invalid = False
        for calc in self._calculations:
            try:
                energy = calc.energy.total.value
            except (IndexError, AttributeError):
                invalid = True
                break
            if energy is None:
                invalid = True
                break
            energies.append(energy.magnitude)
        if invalid:
            logger.warning('Energy not reported for an calculation that is part of a geometry optimization')
        if energies:
            self.results.energies = energies

        if not self.results.final_energy_difference:
            self.results.final_energy_difference = resolve_difference(energies)

        if not self.results.final_force_maximum:
            if len(self._calculations) > 0:
                if self._calculations[-1].forces is not None and self._calculations[-1].forces.total is not None:
                    forces = self._calculations[-1].forces.total.value
                    if forces is not None:
                        max_force = np.max(np.linalg.norm(forces.magnitude, axis=1))
                        self.results.final_force_maximum = max_force * forces.units

        if not self.results.final_displacement_maximum:
            def get_atoms(index):
                system = self._systems[index]
                atoms = Atoms(
                    positions=system.atoms.positions.magnitude,
                    cell=system.atoms.lattice_vectors.magnitude,
                    pbc=system.atoms.periodic)
                atoms.wrap()
                return atoms

            n_systems = len(self._systems)
            a_pos = get_atoms(n_systems - 1).get_positions()
            for i in range(n_systems - 2, -1, -1):
                b_pos = get_atoms(i).get_positions()
                displacement_maximum = np.max(np.abs(a_pos - b_pos))
                if displacement_maximum > 0:
                    self.results.final_displacement_maximum = displacement_maximum
                    break

        if not self.results.is_converged_geometry:
            # we can have several criteria for convergence: energy, force, displacement
            criteria = []
            try:
                criteria.append(self.results.final_energy_difference <= self.method.convergence_tolerance_energy_difference)
            except Exception:
                pass

            try:
                criteria.append(self.results.final_force_maximum <= self.method.convergence_tolerance_force_maximum)
            except Exception:
                pass

            try:
                criteria.append(self.results.final_displacement_maximum <= self.method.convergence_tolerance_displacement_maximum)
            except Exception:
                pass

            # converged when either criterion is met
            if criteria:
                self.results.is_converged_geometry = True in criteria

        if not self.results.calculation_result_ref and self._calculations:
            self.results.calculation_result_ref = self._calculations[-1]


class ThermostatParameters(MSection):
    '''
    Section containing the parameters pertaining to the thermostat for a molecular dynamics run.
    '''

    m_def = Section(validate=False)

    thermostat_type = Quantity(
        type=MEnum('andersen', 'berendsen', 'brownian', 'langevin_goga', 'langevin_schneider', 'nose_hoover', 'velocity_rescaling',
                   'velocity_rescaling_langevin'),
        shape=[],
        description='''
        The name of the thermostat used for temperature control. If skipped or an empty string is used, it
        means no thermostat was applied.

        Allowed values are:

        | Thermostat Name        | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `""`                   | No thermostat               |

        | `"andersen"`           | H.C. Andersen, [J. Chem. Phys.
        **72**, 2384 (1980)](https://doi.org/10.1063/1.439486) |

        | `"berendsen"`          | H. J. C. Berendsen, J. P. M. Postma,
        W. F. van Gunsteren, A. DiNola, and J. R. Haak, [J. Chem. Phys.
        **81**, 3684 (1984)](https://doi.org/10.1063/1.448118) |

        | `"brownian"`           | Brownian Dynamics |

        | `"langevin_goga"`           | N. Goga, A. J. Rzepiela, A. H. de Vries,
        S. J. Marrink, and H. J. C. Berendsen, [J. Chem. Theory Comput. **8**, 3637 (2012)]
        (https://doi.org/10.1021/ct3000876) |

        | `"langevin_schneider"`           | T. Schneider and E. Stoll,
        [Phys. Rev. B **17**, 1302](https://doi.org/10.1103/PhysRevB.17.1302) |

        | `"nose_hoover"`        | S. Nosé, [Mol. Phys. **52**, 255 (1984)]
        (https://doi.org/10.1080/00268978400101201); W.G. Hoover, [Phys. Rev. A
        **31**, 1695 (1985) |

        | `"velocity_rescaling"` | G. Bussi, D. Donadio, and M. Parrinello,
        [J. Chem. Phys. **126**, 014101 (2007)](https://doi.org/10.1063/1.2408420) |

        | `"velocity_rescaling_langevin"` | G. Bussi and M. Parrinello,
        [Phys. Rev. E **75**, 056707 (2007)](https://doi.org/10.1103/PhysRevE.75.056707) |
        ''')

    reference_temperature = Quantity(
        type=np.float64,
        shape=[],
        unit='kelvin',
        description='''
        The target temperature for the simulation.
        ''')

    coupling_constant = Quantity(
        type=np.float64,
        shape=[],
        unit='s',
        description='''
        The time constant for temperature coupling. Need to describe what this means for the various
        thermostat options...
        ''')

    effective_mass = Quantity(
        type=np.float64,
        shape=[],
        unit='kilogram',
        description='''
        The effective or fictitious mass of the temperature resevoir.
        ''')


class BarostatParameters(MSection):
    '''
    Section containing the parameters pertaining to the barostat for a molecular dynamics run.
    '''

    m_def = Section(validate=False)

    barostat_type = Quantity(
        type=MEnum('berendsen', 'martyna_tuckerman_tobias_klein', 'nose_hoover', 'parrinello_rahman', 'stochastic_cell_rescaling'),
        shape=[],
        description='''
        The name of the barostat used for temperature control. If skipped or an empty string is used, it
        means no barostat was applied.

        Allowed values are:

        | Barostat Name          | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `""`                   | No thermostat               |

        | `"berendsen"`          | H. J. C. Berendsen, J. P. M. Postma,
        W. F. van Gunsteren, A. DiNola, and J. R. Haak, [J. Chem. Phys.
        **81**, 3684 (1984)](https://doi.org/10.1063/1.448118) |

        | `"martyna_tuckerman_tobias_klein"` | G.J. Martyna, M.E. Tuckerman, D.J. Tobias, and M.L. Klein,
        [Mol. Phys. **87**, 1117 (1996)](https://doi.org/10.1080/00268979600100761);
        M.E. Tuckerman, J. Alejandre, R. López-Rendón, A.L. Jochim, and G.J. Martyna,
        [J. Phys. A. **59**, 5629 (2006)](https://doi.org/10.1088/0305-4470/39/19/S18)|

        | `"nose_hoover"`        | S. Nosé, [Mol. Phys. **52**, 255 (1984)]
        (https://doi.org/10.1080/00268978400101201); W.G. Hoover, [Phys. Rev. A
        **31**, 1695 (1985) |

        | `"parrinello_rahman"`        | M. Parrinello and A. Rahman,
        [J. Appl. Phys. **52**, 7182 (1981)](https://doi.org/10.1063/1.328693);
        S. Nosé and M.L. Klein, [Mol. Phys. **50**, 1055 (1983) |

        | `"stochastic_cell_rescaling"` | M. Bernetti and G. Bussi,
        [J. Chem. Phys. **153**, 114107 (2020)](https://doi.org/10.1063/1.2408420) |
        ''')

    coupling_type = Quantity(
        type=MEnum('isotropic', 'semi_isotropic', 'anisotropic'),
        shape=[],
        description='''
        Describes the symmetry of pressure coupling. Specifics can be inferred from the `coupling constant`

        | Type          | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `isotropic`          | Identical coupling in all directions. |

        | `semi_isotropic` | Identical coupling in 2 directions. |

        | `anisotropic`        | General case. |
        ''')

    reference_pressure = Quantity(
        type=np.float64,
        shape=[3, 3],
        unit='pascal',
        description='''
        The target pressure for the simulation, stored in a 3x3 matrix, indicating the values for individual directions
        along the diagonal, and coupling between directions on the off-diagonal.
        ''')

    coupling_constant = Quantity(
        type=np.float64,
        shape=[3, 3],
        unit='s',
        description='''
        The time constants for pressure coupling, stored in a 3x3 matrix, indicating the values for individual directions
        along the diagonal, and coupling between directions on the off-diagonal. 0 values along the off-diagonal
        indicate no-coupling between these directions.
        ''')

    compressibility = Quantity(
        type=np.float64,
        shape=[3, 3],
        unit='1 / pascal',
        description='''
        An estimate of the system's compressibility, used for box rescaling, stored in a 3x3 matrix indicating the values for individual directions
        along the diagonal, and coupling between directions on the off-diagonal. If None, it may indicate that these values
        are incorporated into the coupling_constant, or simply that the software used uses a fixed value that is not available in
        the input/output files.
        ''')


class MolecularDynamicsMethod(SimulationWorkflowMethod):

    thermodynamic_ensemble = Quantity(
        type=MEnum('NVE', 'NVT', 'NPT', 'NPH'),
        shape=[],
        description='''
        The type of thermodynamic ensemble that was simulated.

        Allowed values are:

        | Thermodynamic Ensemble          | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `"NVE"`           | Constant number of particles, volume, and energy |

        | `"NVT"`           | Constant number of particles, volume, and temperature |

        | `"NPT"`           | Constant number of particles, pressure, and temperature |

        | `"NPH"`           | Constant number of particles, pressure, and enthalpy |
        ''')

    integrator_type = Quantity(
        type=MEnum(
            'brownian', 'conjugant_gradient', 'langevin_goga',
            'langevin_schneider', 'leap_frog', 'rRESPA_multitimescale', 'velocity_verlet'
        ),
        shape=[],
        description='''
        Name of the integrator.

        Allowed values are:

        | Integrator Name          | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `"langevin_goga"`           | N. Goga, A. J. Rzepiela, A. H. de Vries,
        S. J. Marrink, and H. J. C. Berendsen, [J. Chem. Theory Comput. **8**, 3637 (2012)]
        (https://doi.org/10.1021/ct3000876) |

        | `"langevin_schneider"`           | T. Schneider and E. Stoll,
        [Phys. Rev. B **17**, 1302](https://doi.org/10.1103/PhysRevB.17.1302) |

        | `"leap_frog"`          | R.W. Hockney, S.P. Goel, and J. Eastwood,
        [J. Comp. Phys. **14**, 148 (1974)](https://doi.org/10.1016/0021-9991(74)90010-2) |

        | `"velocity_verlet"` | W.C. Swope, H.C. Andersen, P.H. Berens, and K.R. Wilson,
        [J. Chem. Phys. **76**, 637 (1982)](https://doi.org/10.1063/1.442716) |

        | `"rRESPA_multitimescale"` | M. Tuckerman, B. J. Berne, and G. J. Martyna
        [J. Chem. Phys. **97**, 1990 (1992)](https://doi.org/10.1063/1.463137) |
        ''')

    integration_timestep = Quantity(
        type=np.float64,
        shape=[],
        unit='s',
        description='''
        The timestep at which the numerical integration is performed.
        ''')

    n_steps = Quantity(
        type=int,
        shape=[],
        description='''
        Number of timesteps performed.
        ''')

    coordinate_save_frequency = Quantity(
        type=int,
        shape=[],
        description='''
        The number of timesteps between saving the coordinates.
        ''')

    velocity_save_frequency = Quantity(
        type=int,
        shape=[],
        description='''
        The number of timesteps between saving the velocities.
        ''')

    force_save_frequency = Quantity(
        type=int,
        shape=[],
        description='''
        The number of timesteps between saving the forces.
        ''')

    thermodynamics_save_frequency = Quantity(
        type=int,
        shape=[],
        description='''
        The number of timesteps between saving the thermodynamic quantities.
        ''')

    thermostat_parameters = SubSection(sub_section=ThermostatParameters.m_def, repeats=False)

    barostat_parameters = SubSection(sub_section=BarostatParameters.m_def, repeats=False)


class EnsemblePropertyValues(MSection):
    '''
    Generic section containing information regarding the values of an ensemble property.
    '''

    m_def = Section(validate=False)

    label = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the atoms or molecule types involved in determining the property.
        ''')

    n_bins = Quantity(
        type=int,
        shape=[],
        description='''
        Number of bins.
        ''')

    frame_start = Quantity(
        type=int,
        shape=[],
        description='''
        Trajectory frame number where the ensemble averaging starts.
        ''')

    frame_end = Quantity(
        type=int,
        shape=[],
        description='''
        Trajectory frame number where the ensemble averaging ends.
        ''')


class RadialDistributionFunctionValues(EnsemblePropertyValues):
    '''
    Section containing information regarding the values of
    radial distribution functions (rdfs).
    '''

    m_def = Section(validate=False)

    bins = Quantity(
        type=np.float64,
        shape=['n_bins'],
        unit='m',
        description='''
        Distances along which the rdf was calculated.
        ''')

    value = Quantity(
        type=np.float64,
        shape=['n_bins'],
        description='''
        Values of the property.
        ''')


class EnsembleProperty(ArchiveSection):
    '''
    Generic section containing information about a calculation of any static observable
    from a trajectory (i.e., from an ensemble average).
    '''

    m_def = Section(validate=False)

    type = Quantity(
        type=MEnum('molecular', 'atomic'),
        shape=[],
        description='''
        Describes if the observable is calculated at the molecular or atomic level.
        ''')

    n_smooth = Quantity(
        type=int,
        shape=[],
        description='''
        Number of bins over which the running average was computed for
        the observable `values'.
        ''')

    error_type = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the type of error reported for this observable.
        ''')

    n_variables = Quantity(
        type=int,
        shape=[],
        description='''
        Number of variables along which the property is determined.
        ''')

    variables_name = Quantity(
        type=str,
        shape=['n_variables'],
        description='''
        Name/description of the independent variables along which the observable is defined.
        ''')


class RadialDistributionFunction(EnsembleProperty):
    '''
    Section containing information about the calculation of
    radial distribution functions (rdfs).
    '''

    m_def = Section(validate=False)

    _rdf_results = None

    radial_distribution_function_values = SubSection(sub_section=RadialDistributionFunctionValues.m_def, repeats=True)

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if self._rdf_results:
            self.type = self._rdf_results.get('type')
            self.n_smooth = self._rdf_results.get('n_smooth')
            self.n_prune = self._rdf_results.get('n_prune')
            self.n_variables = 1
            self.variables_name = ['distance']
            for i_pair, pair_type in enumerate(self._rdf_results.get('types', [])):
                sec_rdf_values = self.m_create(RadialDistributionFunctionValues)
                sec_rdf_values.label = str(pair_type)
                sec_rdf_values.n_bins = len(self._rdf_results.get('bins', [[]] * i_pair)[i_pair])
                sec_rdf_values.bins = self._rdf_results.get('bins', [[]] * i_pair)[i_pair]
                sec_rdf_values.value = self._rdf_results.get('value', [[]] * i_pair)[i_pair]
                sec_rdf_values.frame_start = self._rdf_results.get('frame_start', [[]] * i_pair)[i_pair]
                sec_rdf_values.frame_end = self._rdf_results.get('frame_end', [[]] * i_pair)[i_pair]


class TrajectoryProperty(ArchiveSection):
    '''
    Generic section containing information about a calculation of any observable
    defined and stored at each individual frame of a trajectory.
    '''

    m_def = Section(validate=False)

    type = Quantity(
        type=MEnum('molecular', 'atomic'),
        shape=[],
        description='''
        Describes if the observable is calculated at the molecular or atomic level.
        ''')

    error_type = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the type of error reported for this observable.
        ''')

    label = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the atoms or molecule types involved in determining the property.
        ''')

    n_frames = Quantity(
        type=int,
        shape=[],
        description='''
        Number of frames for which the observable is stored.
        ''')

    frames = Quantity(
        type=np.int32,
        shape=['n_frames'],
        description='''
        Frames for which the observable is stored.
        ''')

    times = Quantity(
        type=np.float64,
        shape=['n_frames'],
        unit='s',
        description='''
        Times for which the observable is stored.
        ''')

    value = Quantity(
        type=np.float64,
        shape=['n_frames'],
        description='''
        Values of the property.
        ''')

    errors = Quantity(
        type=np.float64,
        shape=['*'],
        description='''
        Error associated with the determination of the property.
        ''')


class RadiusOfGyration(TrajectoryProperty):
    '''
    Section containing information about the calculation of
    radius of gyration (Rg).
    '''

    m_def = Section(validate=False)

    _rg_results = None

    atomsgroup_ref = Quantity(
        type=Reference(AtomsGroup.m_def),
        shape=[1],
        description='''
        References to the atoms_group section containing the molecule for which Rg was calculated.
        ''')

    value = Quantity(
        type=np.float64,
        shape=['n_frames'],
        unit='m',
        description='''
        Values of the property.
        ''')

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if self._rg_results:
            self.type = self._rg_results.get('type')
            self.label = self._rg_results.get('label')
            self.atomsgroup_ref = self._rg_results.get('atomsgroup_ref')
            self.n_frames = self._rg_results.get('n_frames')
            self.times = self._rg_results.get('times')
            self.value = self._rg_results.get('value')


class DiffusionConstantValues(MSection):
    '''
    Section containing information regarding the diffusion constants.
    '''

    m_def = Section(validate=False)

    value = Quantity(
        type=np.float64,
        shape=[],
        unit='m^2/s',
        description='''
        Values of the diffusion constants.
        ''')

    error_type = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the type of error reported for this observable.
        ''')

    errors = Quantity(
        type=np.float64,
        shape=['*'],
        description='''
        Error associated with the determination of the diffusion constant.
        ''')


class CorrelationFunctionValues(MSection):
    '''
    Generic section containing information regarding the values of a correlation function.
    '''

    m_def = Section(validate=False)

    label = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the atoms or molecule types involved in determining the property.
        ''')

    n_times = Quantity(
        type=int,
        shape=[],
        description='''
        Number of times windows for the calculation of the correlation function.
        ''')


class MeanSquaredDisplacementValues(CorrelationFunctionValues):
    '''
    Section containing information regarding the values of a mean squared displacements (msds).
    '''

    m_def = Section(validate=False)

    times = Quantity(
        type=np.float64,
        shape=['n_times'],
        unit='s',
        description='''
        Time windows used for the calculation of the msds.
        ''')

    value = Quantity(
        type=np.float64,
        shape=['n_times'],
        unit='m^2',
        description='''
        Msd values.
        ''')

    errors = Quantity(
        type=np.float64,
        shape=['*'],
        description='''
        Error associated with the determination of the msds.
        ''')

    diffusion_constant = SubSection(sub_section=DiffusionConstantValues.m_def, repeats=False)


class CorrelationFunction(ArchiveSection):
    '''
    Generic section containing information about a calculation of any time correlation
    function from a trajectory.
    '''

    m_def = Section(validate=False)

    type = Quantity(
        type=MEnum('molecular', 'atomic'),
        shape=[],
        description='''
        Describes if the correlation function is calculated at the molecular or atomic level.
        ''')

    direction = Quantity(
        type=MEnum('x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'),
        shape=[],
        description='''
        Describes the direction in which the correlation function was calculated.
        ''')

    error_type = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the type of error reported for this correlation function.
        ''')


class MeanSquaredDisplacement(CorrelationFunction):
    '''
    Section containing information about a calculation of any mean squared displacements (msds).
    '''

    m_def = Section(validate=False)

    _msd_results = None

    mean_squared_displacement_values = SubSection(sub_section=MeanSquaredDisplacementValues.m_def, repeats=True)

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if not self._msd_results:
            return

        self.type = self._msd_results.get('type')
        self.direction = self._msd_results.get('direction')
        for i_type, moltype in enumerate(self._msd_results.get('types', [])):
            sec_msd_values = self.m_create(MeanSquaredDisplacementValues)
            sec_msd_values.label = str(moltype)
            sec_msd_values.n_times = len(self._msd_results.get('times', [[]] * i_type)[i_type])
            sec_msd_values.times = self._msd_results['times'][i_type] if self._msd_results.get(
                'times') is not None else []
            sec_msd_values.value = self._msd_results['value'][i_type] if self._msd_results.get(
                'value') is not None else []
            sec_diffusion = sec_msd_values.m_create(DiffusionConstantValues)
            sec_diffusion.value = self._msd_results['diffusion_constant'][i_type] if self._msd_results.get(
                'diffusion_constant') is not None else []
            sec_diffusion.error_type = 'Pearson correlation coefficient'
            sec_diffusion.errors = self._msd_results['error_diffusion_constant'][i_type] if self._msd_results.get(
                'error_diffusion_constant') is not None else []


class MolecularDynamicsResults(ThermodynamicsResults):

    finished_normally = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if calculation terminated normally.
        ''')

    n_steps = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of trajectory steps''')

    trajectory = Quantity(
        type=Reference(System),
        shape=['n_steps'],
        description='''
        Reference to the system of each step in the trajectory.
        ''')

    radial_distribution_functions = SubSection(sub_section=RadialDistributionFunction, repeats=True)

    radius_of_gyration = SubSection(sub_section=RadiusOfGyration, repeats=True)

    mean_squared_displacements = SubSection(sub_section=MeanSquaredDisplacement, repeats=True)

    def normalize(self, archive, logger):

        super().normalize(archive, logger)

        universe = archive_to_universe(archive)
        if universe is None:
            return

        # calculate molecular radial distribution functions
        if not self.radial_distribution_functions:
            n_traj_split = 10  # number of intervals to split trajectory into for averaging
            interval_indices = []  # 2D array specifying the groups of the n_traj_split intervals to be averaged
            # first 20% of trajectory
            interval_indices.append(np.arange(int(n_traj_split * 0.20)))
            # last 80% of trajectory
            interval_indices.append(np.arange(n_traj_split)[len(interval_indices[0]):])
            # last 60% of trajectory
            interval_indices.append(np.arange(n_traj_split)[len(interval_indices[0]) * 2:])
            # last 40% of trajectory
            interval_indices.append(np.arange(n_traj_split)[len(interval_indices[0]) * 3:])

            n_prune = int(universe.trajectory.n_frames / len(archive.run[-1].system))
            rdf_results = calc_molecular_rdf(universe, n_traj_split=n_traj_split,
                                             n_prune=n_prune, interval_indices=interval_indices)
            if rdf_results:
                sec_rdfs = RadialDistributionFunction()
                sec_rdfs._rdf_results = rdf_results
                self.radial_distribution_functions.append(sec_rdfs)

        # calculate the molecular mean squared displacements
        if not self.mean_squared_displacements:
            msd_results = calc_molecular_mean_squared_displacements(universe)
            if msd_results:
                sec_msds = MeanSquaredDisplacement()
                sec_msds._msd_results = msd_results
                self.mean_squared_displacements.append(sec_msds)

        # calculate radius of gyration for polymers
        try:
            sec_system = archive.run[-1].system[0]
            sec_calc = archive.run[-1].calculation
            sec_calc = sec_calc if sec_calc is not None else []
        except Exception:
            return

        flag_rgs = False
        for i_calc, calc in enumerate(sec_calc):
            if calc.radius_of_gyration:
                flag_rgs = True
                break  # TODO Should transfer Rg's to workflow results if they are already supplied in calculation

        if not flag_rgs:
            flag_warned = False
            sec_rgs_calc = None
            system_topology = sec_system.get('atoms_group')
            rg_results = calc_molecular_radius_of_gyration(universe, system_topology)
            for rg in rg_results:
                sec_rgs = RadiusOfGyration()
                sec_rgs._rg_results = rg
                self.radius_of_gyration.append(sec_rgs)

                # disperse results into calculation section
                n_frames = rg.get('n_frames')
                if len(sec_calc) == 0:
                    sec_calc = [self.run.m_create(Calculation) for _ in range(n_frames)]
                elif n_frames != len(sec_calc):
                    if not flag_warned:
                        self.logger.warning(
                            'Unexpected mismatch in number of calculations and number of'
                            'trajectory frames. Not storing Rg values under calculation.')
                        flag_warned = True
                    # TODO sync calculation and system sections to be able to store the Rgs under calculation
                    continue
                for i_calc, calc in enumerate(sec_calc):
                    sec_rgs_calc = calc.radius_of_gyration
                    if not sec_rgs_calc:
                        sec_rgs_calc = calc.m_create(RadiusOfGyrationCalculation)
                        sec_rgs_calc.kind = rg.get('type')
                    else:
                        sec_rgs_calc = sec_rgs_calc[0]
                    sec_rg_values = sec_rgs_calc.m_create(RadiusOfGyrationValuesCalculation)
                    sec_rg_values.atomsgroup_ref = rg.get('atomsgroup_ref')
                    sec_rg_values.label = rg.get('label')
                    sec_rg_values.value = rg.get('value')[i_calc]


class MolecularDynamics(SerialSimulation):

    method = SubSection(sub_section=MolecularDynamicsMethod)

    results = SubSection(sub_section=MolecularDynamicsResults)

    def normalize(self, archive: EntryArchive, logger):

        super().normalize(archive, logger)

        if not self.method:
            self.method = MolecularDynamicsMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = MolecularDynamicsResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))

        if self.results.trajectory is None and self._systems:
            self.results.trajectory = self._systems

        def set_thermo_property(name):
            values = []
            quantity = None
            for calc in self._calculations:
                try:
                    quantity = getattr(calc, name)
                    if quantity is not None:
                        values.append(quantity.magnitude if hasattr(quantity, 'magnitude') else quantity)
                except Exception:
                    pass
            if len(values) == 0:
                return

            unit = quantity.units if hasattr(quantity, 'units') else 1.0
            setattr(self.results, name, np.array(values) * unit)

        if self.results.temperature is None:
            set_thermo_property('temperature')

        if self.results.pressure is None:
            set_thermo_property('pressure')

        if not self.results.calculation_result_ref and self._calculations:
            self.results.calculation_result_ref = self._calculations[-1]


class PhononMethod(SimulationWorkflowMethod):

    force_calculator = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the program used to calculate the forces.
        ''')

    mesh_density = Quantity(
        type=np.float64,
        shape=[],
        unit='1 / meter ** 3',
        description='''
        Density of the k-mesh for sampling.
        ''')

    random_displacements = Quantity(
        type=bool,
        shape=[],
        description='''
        Identifies if displacements are made randomly.
        ''')

    with_non_analytic_correction = Quantity(
        type=bool,
        shape=[],
        description='''
        Identifies if non-analytical term corrections are applied to dynamical matrix.
        ''')

    with_grueneisen_parameters = Quantity(
        type=bool,
        shape=[],
        description='''
        Identifies if Grueneisen parameters are calculated.
        ''')


class PhononResults(ThermodynamicsResults):

    n_imaginary_frequencies = Quantity(
        type=int,
        shape=[],
        description='''
        Number of modes with imaginary frequencies.
        ''')

    n_bands = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of phonon bands.
        ''')

    n_qpoints = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of q points for which phonon properties are evaluated.
        ''')

    qpoints = Quantity(
        type=np.float64,
        shape=['n_qpoints', 3],
        description='''
        Value of the qpoints.
        ''')

    group_velocity = Quantity(
        type=np.float64,
        shape=['n_qpoints', 'n_bands', 3],
        unit='meter / second',
        description='''
        Calculated value of the group velocity at each qpoint.
        ''')

    n_displacements = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of independent displacements.
        ''')

    n_atoms = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of atoms in the simulation cell.
        ''')

    displacements = Quantity(
        type=np.float64,
        shape=['n_displacements', 'n_atoms', 3],
        unit='meter',
        description='''
        Value of the displacements applied to each atom in the simulation cell.
        ''')

    dos = Quantity(
        type=Reference(Dos),
        shape=['n_data'],
        description='''
        Reference to the electronic density of states data.
        ''')

    band_structure = Quantity(
        type=Reference(BandStructure),
        shape=['n_data'],
        description='''
        Reference to the electronic band structure data.
        ''')


class Phonon(ParallelSimulation):

    method = SubSection(sub_section=PhononMethod)

    results = SubSection(sub_section=PhononResults)

    def normalize(self, archive: EntryArchive, logger):

        super().normalize(archive, logger)

        if not self._calculations or not self._calculations[0].band_structure_phonon:
            return

        if not self.method:
            self.method = PhononMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = PhononResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))

        last_calc = self._calculations[-1]

        if not self.results.n_imaginary_frequencies:
            n_imaginary = 0
            for band_segment in last_calc.band_structure_phonon[-1].segment:
                freq = band_segment.energies.magnitude
                n_imaginary += np.count_nonzero(np.array(freq) < 0)
            self.results.n_imaginary_frequencies = n_imaginary

        if not self.results.calculation_result_ref:
            self.results.calculation_result_ref = last_calc

        if not self.results.dos:
            self.results.dos = last_calc.dos_phonon

        if not self.results.band_structure:
            self.results.band_structure = last_calc.band_structure_phonon


class StrainDiagrams(MSection):
    '''
    Section containing the information regarding the elastic strains.
    '''

    m_def = Section(
        validate=False)

    type = Quantity(
        type=str,
        shape=[],
        description='''
        Kind of strain diagram. Possible values are: energy; cross-validation (cross-
        validation error); d2E (second derivative of the energy wrt the strain)
        ''')

    n_eta = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of strain values used in the strain diagram
        ''')

    n_deformations = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of deformations.
        ''')

    value = Quantity(
        type=np.float64,
        shape=['n_deformations', 'n_eta'],
        description='''
        Values of the energy(units:J)/d2E(units:Pa)/cross-validation (depending on the
        value of strain_diagram_type)
        ''')

    eta = Quantity(
        type=np.float64,
        shape=['n_deformations', 'n_eta'],
        description='''
        eta values used the strain diagrams
        ''')

    stress_voigt_component = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Voigt component corresponding to the strain diagram
        ''')

    polynomial_fit_order = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Order of the polynomial fit
        ''')


class ElasticMethod(SimulationWorkflowMethod):

    energy_stress_calculator = Quantity(
        type=str,
        shape=[],
        description='''
        Name of program used to calculate energy or stress.
        ''')

    calculation_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to calculate elastic constants, can either be energy or stress.
        ''')

    elastic_constants_order = Quantity(
        type=int,
        shape=[],
        description='''
        Order of the calculated elastic constants.
        ''')

    fitting_error_maximum = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Maximum error in polynomial fit.
        ''')

    strain_maximum = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Maximum strain applied to crystal.
        ''')


class ElasticResults(ThermodynamicsResults):

    n_deformations = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of deformed structures used to calculate the elastic constants. This is
        determined by the symmetry of the crystal.
        ''')

    deformation_types = Quantity(
        type=np.str_,
        shape=['n_deformations', 6],
        description='''
        deformation types
        ''')

    n_strains = Quantity(
        type=np.int32,
        shape=[],
        description='''
        number of equally spaced strains applied to each deformed structure, which are
        generated between the maximum negative strain and the maximum positive one.
        ''')

    is_mechanically_stable = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if structure is mechanically stable from the calculated values of the
        elastic constants.
        ''')

    elastic_constants_notation_matrix_second_order = Quantity(
        type=np.str_,
        shape=[6, 6],
        description='''
        Symmetry of the second-order elastic constant matrix in Voigt notation
        ''')

    elastic_constants_matrix_second_order = Quantity(
        type=np.float64,
        shape=[6, 6],
        unit='pascal',
        description='''
        2nd order elastic constant (stiffness) matrix in pascals
        ''')

    elastic_constants_matrix_third_order = Quantity(
        type=np.float64,
        shape=[6, 6, 6],
        unit='pascal',
        description='''
        3rd order elastic constant (stiffness) matrix in pascals
        ''')

    compliance_matrix_second_order = Quantity(
        type=np.float64,
        shape=[6, 6],
        unit='1 / pascal',
        description='''
        Elastic compliance matrix in 1/GPa
        ''')

    elastic_constants_gradient_matrix_second_order = Quantity(
        type=np.float64,
        shape=[18, 18],
        unit='newton',
        description='''
        gradient of the 2nd order elastic constant (stiffness) matrix in newton
        ''')

    bulk_modulus_voigt = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Voigt bulk modulus
        ''')

    shear_modulus_voigt = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Voigt shear modulus
        ''')

    bulk_modulus_reuss = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Reuss bulk modulus
        ''')

    shear_modulus_reuss = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Reuss shear modulus
        ''')

    bulk_modulus_hill = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Hill bulk modulus
        ''')

    shear_modulus_hill = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Hill shear modulus
        ''')

    young_modulus_voigt = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Voigt Young modulus
        ''')

    poisson_ratio_voigt = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Voigt Poisson ratio
        ''')

    young_modulus_reuss = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Reuss Young modulus
        ''')

    poisson_ratio_reuss = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Reuss Poisson ratio
        ''')

    young_modulus_hill = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Hill Young modulus
        ''')

    poisson_ratio_hill = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Hill Poisson ratio
        ''')

    elastic_anisotropy = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Elastic anisotropy
        ''')

    pugh_ratio_hill = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Pugh ratio defined as the ratio between the shear modulus and bulk modulus
        ''')

    debye_temperature = Quantity(
        type=np.float64,
        shape=[],
        unit='kelvin',
        description='''
        Debye temperature
        ''')

    speed_sound_transverse = Quantity(
        type=np.float64,
        shape=[],
        unit='meter / second',
        description='''
        Speed of sound along the transverse direction
        ''')

    speed_sound_longitudinal = Quantity(
        type=np.float64,
        shape=[],
        unit='meter / second',
        description='''
        Speed of sound along the longitudinal direction
        ''')

    speed_sound_average = Quantity(
        type=np.float64,
        shape=[],
        unit='meter / second',
        description='''
        Average speed of sound
        ''')

    eigenvalues_elastic = Quantity(
        type=np.float64,
        shape=[6],
        unit='pascal',
        description='''
        Eigenvalues of the stiffness matrix
        ''')

    strain_diagrams = SubSection(
        sub_section=StrainDiagrams.m_def,
        repeats=True)


class Elastic(ParallelSimulation):

    method = SubSection(sub_section=ElasticMethod)

    results = SubSection(sub_section=ElasticResults)

    def _resolve_mechanical_stability(self):
        spacegroup, c = None, None
        try:
            spacegroup = self._systems[-1].symmetry[-1].space_group_number
            c = self.results.elastic_constants_matrix_second_order
        except Exception:
            return False

        if c is None or spacegroup is None:
            return False

        # see Phys. Rev B 90, 224104 (2014)
        res = False
        if spacegroup <= 2:  # Triclinic
            res = np.count_nonzero(c < 0)
        elif spacegroup <= 15:  # Monoclinic
            res = np.count_nonzero(c < 0)
        elif spacegroup <= 74:  # Orthorhombic
            res =\
                c[0][0] > 0 and c[0][0] * c[1][1] > c[0][1] ** 2 and\
                c[0][0] * c[1][1] * c[2][2] + 2 * c[0][1] * c[0][2] * c[1][2] -\
                c[0][0] * c[1][2] ** 2 - c[1][1] * c[0][2] ** 2 - c[2][2] * c[0][1] ** 2 > 0 and\
                c[3][3] > 0 and c[4][4] > 0 and c[5][5] > 0
        elif spacegroup <= 88:  # Tetragonal II
            res =\
                c[0][0] > abs(c[0][1]) and\
                2 * c[0][2] ** 2 < c[2][2] * (c[0][0] + c[0][1])
        elif spacegroup <= 142:  # Tetragonal I
            res =\
                c[0][0] > abs(c[0][1]) and\
                2 * c[0][2] ** 2 < c[2][2] * (c[0][0] + c[0][1]) and\
                c[3][3] > 0 and c[5][5] > 0
        elif spacegroup <= 148:  # rhombohedral II
            res =\
                c[0][0] > abs(c[0][1]) and c[3][3] > 0 and\
                c[0][2] ** 2 < (0.5 * c[2][2] * (c[0][0] + c[0][1])) and\
                c[0][3] ** 2 + c[0][4] ** 2 < 0.5 * c[3][3] * (c[0][0] - c[0][1])
        elif spacegroup <= 167:  # rhombohedral I
            res =\
                c[0][0] > abs(c[0][1]) and c[3][3] > 0 and\
                c[0][2] ** 2 < 0.5 * c[2][2] * (c[0][0] + c[0][1]) and\
                c[0][3] ** 2 < 0.5 * c[3][3] * (c[0][0] - c[0][1])
        elif spacegroup <= 194:  # hexagonal I
            res =\
                c[0][0] > abs(c[0][1]) and\
                2 * c[0][2] ** 2 < c[2][2] * (c[0][0] + c[0][1]) and\
                c[3][3] > 0 and c[5][5] > 0
        else:  # cubic
            res = c[0][0] - c[0][1] > 0 and c[0][0] + 2 * c[0][1] > 0 and c[3][3] > 0

        return res

    def _get_maximum_fit_error(self):
        max_error = 0.0
        if len(self._calculations) == 0:
            return max_error

        for diagram in self.results.strain_diagrams:
            if diagram.type == 'cross-validation':
                error = np.amax(diagram.value)
                max_error = error if error > max_error else max_error

        return max_error

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = ElasticMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = ElasticResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))

        if self.results.is_mechanically_stable is None:
            self.results.is_mechanically_stable = bool(self._resolve_mechanical_stability())

        if self.method.fitting_error_maximum is None:
            self.method.fitting_error_maximum = self._get_maximum_fit_error()

        if not self.results.calculation_result_ref and self._calculations:
            self.results.calculation_result_ref = self._calculations[-1]


class ThermodynamicsMethod(SimulationWorkflowMethod):

    pass


class Thermodynamics(SerialSimulation):

    method = SubSection(sub_section=ThermodynamicsMethod)

    results = SubSection(sub_section=ThermodynamicsResults)

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = ThermodynamicsMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = ThermodynamicsResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))


class GWResults(SimulationWorkflowResults):

    dos_dft = Quantity(
        type=Reference(Dos),
        description='''
        DFT density of states
        ''')

    dos_gw = Quantity(
        type=Reference(Dos),
        description='''
        GW density of states
        ''')

    band_structure_dft = Quantity(
        type=Reference(BandStructure),
        description='''
        DFT density of states
        ''')

    band_structure_gw = Quantity(
        type=Reference(BandStructure),
        description='''
        DFT density of states
        ''')


class GWMethod(SimulationWorkflowMethod):

    gw_method_ref = Quantity(
        type=Reference(GWMethodology),
        description='''
        GW methodology reference.
        ''')

    starting_point = Quantity(
        type=Reference(XCFunctional),
        description='''
        Starting point (XC functional or HF) used.
        ''')

    basis_set = Quantity(
        type=Reference(BasisSet),
        description='''
        Basis set used.
        ''')


class GW(SerialSimulation):

    method = SubSection(sub_section=GWMethod)

    results = SubSection(sub_section=GWResults)

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = GWMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))
            # link method also to first task
            if self.tasks:
                self.tasks[0].inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = GWResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))
            # link results also to last task
            if self.tasks:
                self.tasks[-1].inputs.append(Link(name=_workflow_results_name, section=self.results))


class PhotonPolarizationResults(SimulationWorkflowResults):

    n_polarizations = Quantity(
        type=np.int32,
        description='''
        Number of polarizations for the phonons used for the calculations.
        ''')

    spectrum_polarization = Quantity(
        type=Reference(Spectra),
        shape=['n_polarizations'],
        description='''
        Spectrum for a given polarization of the photon.
        ''')


class PhotonPolarizationMethod(SimulationWorkflowMethod):

    pass


class PhotonPolarization(ParallelSimulation):

    method = SubSection(sub_section=PhotonPolarizationMethod)

    results = SubSection(sub_section=PhotonPolarizationResults)

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = PhotonPolarizationMethod()
            if not self.inputs:
                self.inputs.append(Link(name=_workflow_method_name, section=self.method))
                # link method also to first task
                if self.tasks:
                    self.tasks[0].inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = PhotonPolarizationResults()
            if not self.outputs:
                self.outputs.append(Link(name=_workflow_results_name, section=self.results))
                # link results also to last task
                if self.tasks:
                    self.tasks[-1].inputs.append(Link(name=_workflow_results_name, section=self.results))


class ParticleHoleExcitationsResults(SimulationWorkflowResults):

    dos_dft = Quantity(
        type=Reference(Dos),
        description='''
        DFT density of states
        ''')

    dos_gw = Quantity(
        type=Reference(Dos),
        description='''
        GW density of states
        ''')

    band_structure_dft = Quantity(
        type=Reference(BandStructure),
        description='''
        DFT density of states
        ''')

    band_structure_gw = Quantity(
        type=Reference(BandStructure),
        description='''
        DFT density of states
        ''')

    spectra = SubSection(sub_section=PhotonPolarizationResults, repeats=True)


class ParticleHoleExcitationsMethod(SimulationWorkflowMethod):

    pass


class ParticleHoleExcitations(SerialSimulation):

    method = SubSection(sub_section=ParticleHoleExcitationsMethod)

    results = SubSection(sub_section=ParticleHoleExcitationsResults)

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = ParticleHoleExcitationsMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))
            # link method also to first task
            if self.tasks:
                self.tasks[0].inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = ParticleHoleExcitationsResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))
            # link results also to last task
            if self.tasks:
                self.tasks[-1].inputs.append(Link(name=_workflow_results_name, section=self.results))


class EquationOfStateMethod(SimulationWorkflowMethod):

    energy_calculator = Quantity(
        type=str,
        shape=[],
        description='''
        Name of program used to calculate energy.
        ''')


class EOSFit(MSection):
    '''
    Section containing results of an equation of state fit.
    '''

    m_def = Section(validate=False)

    function_name = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the function used to perform the fitting of the volume-energy data. Value
        can be one of birch_euler, birch_lagrange, birch_murnaghan, mie_gruneisen,
        murnaghan, pack_evans_james, poirier_tarantola, tait, vinet.
        ''')

    fitted_energies = Quantity(
        type=np.float64,
        shape=['n_points'],
        unit='joule',
        description='''
        Array of the fitted energies corresponding to each volume.
        ''')

    bulk_modulus = Quantity(
        type=np.float64,
        shape=[],
        unit='pascal',
        description='''
        Calculated value of the bulk modulus by fitting the volume-energy data.
        ''')

    bulk_modulus_derivative = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Calculated value of the pressure derivative of the bulk modulus.
        ''')

    equilibrium_volume = Quantity(
        type=np.float64,
        shape=[],
        unit='m ** 3',
        description='''
        Calculated value of the equilibrium volume by fitting the volume-energy data.
        ''')

    equilibrium_energy = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Calculated value of the equilibrium energy by fitting the volume-energy data.
        ''')

    rms_error = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Root-mean squared value of the error in the fitting.
        ''')


class EquationOfStateResults(SimulationWorkflowResults):

    n_points = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of volume-energy pairs in data.
        ''')

    volumes = Quantity(
        type=np.float64,
        shape=['n_points'],
        unit='m ** 3',
        description='''
        Array of volumes per atom for which the energies are evaluated.
        ''')

    energies = Quantity(
        type=np.float64,
        shape=['n_points'],
        unit='joule',
        description='''
        Array of energies corresponding to each volume.
        ''')

    eos_fit = SubSection(sub_section=EOSFit.m_def, repeats=True)


class EquationOfState(ParallelSimulation):

    method = SubSection(sub_section=EquationOfStateMethod)

    results = SubSection(sub_section=EquationOfStateResults)

    def normalize(self, archive: EntryArchive, logger):
        super().normalize(archive, logger)

        if not self.method:
            self.method = EquationOfStateMethod()
            self.inputs.append(Link(name=_workflow_method_name, section=self.method))

        if not self.results:
            self.results = EquationOfStateResults()
            self.outputs.append(Link(name=_workflow_results_name, section=self.results))

        if not self._calculations:
            return

        if not self.results.energies:
            try:
                self.results.energies = [calc.energy.total.value.magnitude for calc in self._calculations]
            except Exception:
                pass

        if not self.results.volumes:
            try:
                self.results.volumes = []
                for system in self._systems:
                    cell = system.atoms.simulation_cell
                    if cell:
                        self.results.volumes.append(np.dot(cell[0], np.cross(cell[1], cell[2])))
            except Exception:
                pass
