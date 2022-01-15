/*
 e Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Contains GUI test data for a typical DFT entry with several different
 * properties:
 *
 * - Electronic DOS
 * - Electronic band structure
 */
import { cloneDeep } from 'lodash'

// Indexed data that is shared between metadata and results
const common = {
  domain: 'dft',
  upload_id: 'mock',
  upload_create_time: '2021-03-17T13:47:32.899000',
  last_processing_time: '2021-03-17T15:47:32.899000',
  nomad_version: '0.10.0',
  nomad_commit: 'bf3c06fa',
  entry_id: 'dft_bulk',
  comment: 'Mocked',
  references: ['doi'],
  authors: [{name: 'Lauri Himanen'}],
  datasets: [{dataset_id: 'Mock dataset', dataset_name: 'Mock dataset'}],
  mainfile: 'vasp.xml',
  formula: 'Si2',
  quantities: []
}

const vdwMethod = 'G06'
const relativityMethod = 'scalar_relativistic_atomic_ZORA'
const basisSetName = 'STO-3G'
const basisSetType = 'plane waves'
const programName = 'VASP'
const programVersion = '1'

const workflow = [
  {
    thermodynamics: [
      {
        heat_capacity_c_v: [0, 1],
        vibrational_free_energy_at_constant_volume: [0, 1],
        temperature: [0, 100]
      }
    ],
    type: 'phonon',
    calculation_result_ref: '/run/0/calculation/0'
  },
  {
    equation_of_state: {
      energies: [0, 1],
      volumes: [0, 1],
      eos_fit: [
        {
          function_name: 'murnaghan',
          fitted_energies: [0, 1],
          bulk_modulus: 1
        }
      ]
    },
    type: 'equation_of_state',
    calculation_result_ref: '/run/0/calculation/0'
  },
  {
    elastic: [
      {
        shear_modulus_hill: 1
      }
    ],
    type: 'elastic',
    calculation_result_ref: '/run/0/calculation/0'
  }
]

// Indexed data that is specific to results
const resultsDftBulk = {
  material: {
    material_id: 'Mock material id',
    material_name: 'Silicon',
    structural_type: 'bulk',
    elements: ['Si'],
    n_elements: 1,
    chemical_formula_reduced: 'Si2',
    chemical_formula_hill: 'Si2',
    chemical_formula_anonymous: 'A2',
    chemical_formula_descriptive: 'Si2',
    symmetry: {
      crystal_system: 'cubic',
      bravais_lattice: 'cP',
      structure_name: 'rock salt',
      space_group_symbol: 'Fd-3m',
      space_group_number: 227,
      point_group: '6mm'
    }
  },
  method: {
    method_name: 'DFT',
    simulation: {
      program_name: programName,
      program_version: programVersion,
      dft: {
        basis_set_type: basisSetType,
        basis_set_name: basisSetName,
        van_der_Waals_method: vdwMethod,
        relativity_method: relativityMethod,
        xc_functional_type: 'GGA',
        xc_functional_names: ['GGA_C_PBE', 'GGA_X_PBE']
      }
    }
  },
  properties: {
    available_properties: [
      'dos_electronic',
      'band_structure_electronic',
      'dos_phonon',
      'band_structure_phonon',
      'heat_capacity_constant_volume',
      'energy_free_helmholtz',
      'bulk_modulus',
      'shear_modulus',
      'energy_volume_curve'
    ],
    electronic: {
      dos_electronic: {
        energies: '/run/0/calculation/0/dos_electronic/0/energies',
        total: '/run/0/calculation/0/dos_electronic/0/total',
        band_gap: [{
          energy_highest_occupied: 0
        }]
      },
      band_structure_electronic: {
        segment: ['/run/0/calculation/0/band_structure_electronic/0/segment/0'],
        reciprocal_cell: '/run/0/calculation/0/band_structure_electronic/0/reciprocal_cell',
        band_gap: [{
          energy_highest_occupied: 0,
          value: 1e-19,
          type: 'indirect'
        }]
      }
    },
    structures: {
      structure_original: {
        lattice_parameters: {
          a: 5e-10,
          b: 5e-10,
          c: 5e-10,
          alpha: Math.PI / 2,
          beta: Math.PI / 2,
          gamma: Math.PI / 2
        },
        cell_volume: 125e-30
      },
      structure_conventional: {
        lattice_parameters: {
          a: 5e-10,
          b: 5e-10,
          c: 5e-10,
          alpha: Math.PI / 2,
          beta: Math.PI / 2,
          gamma: Math.PI / 2
        },
        cell_volume: 125e-37
      },
      structure_primitive: {
        lattice_parameters: {
          a: 5e-10,
          b: 5e-10,
          c: 5e-10,
          alpha: Math.PI / 2,
          beta: Math.PI / 2,
          gamma: Math.PI / 2
        },
        cell_volume: 125e-30
      }
    },
    vibrational: {
      dos_phonon: {
        energies: '/run/0/calculation/0/dos_phonon/0/energies',
        total: '/run/0/calculation/0/dos_phonon/0/total'
      },
      band_structure_phonon: {
        segment: ['/run/0/calculation/0/band_structure_phonon/0/segment/0']
      },
      heat_capacity_constant_volume: {
        heat_capacities: '/workflow/0/thermodynamics/0/heat_capacity_c_v',
        temperatures: '/workflow/0/thermodynamics/0/temperature'
      },
      energy_free_helmholtz: {
        energies: '/workflow/0/thermodynamics/0/vibrational_free_energy_at_constant_volume',
        temperatures: '/workflow/0/thermodynamics/0/temperature'
      }
    },
    mechanical: {
      bulk_modulus: [
        {
          type: 'murnaghan',
          value: 1
        }
      ],
      shear_modulus: [
        {
          type: 'voigt_reuss_hill_average',
          value: 1
        }
      ],
      energy_volume_curve: [
        {
          type: 'murhaghan',
          volumes: '/workflow/1/equation_of_state/volumes',
          energies_fit: '/workflow/1/equation_of_state/eos_fit/0/fitted_energies'
        }
      ]
    }
  }
}

// Section run
const run = [{
  program: {
    name: programName,
    version: programVersion
  },
  method: [
    {
      dft: {
        xc_functional: {
          correlation: [{name: 'GGA_C_PBE'}],
          exchange: [{name: 'GGA_X_PBE'}]
        }
      },
      electronic: {
        van_der_Waals_method: vdwMethod,
        relativity_method: relativityMethod,
        method: 'DFT'
      },
      basis_set: [
        {
          name: basisSetName,
          type: basisSetType
        }
      ]
    }
  ],
  calculation: [
    {
      dos_electronic: [
        {
          energies: [0, 1e-19],
          total: [
            {
              value: [0, 1e18],
              normalization_factor: 1e-19,
              spin: 0
            }
          ],
          band_gap: [{
            energy_highest_occupied: 0,
            index: 0
          }]
        }
      ],
      dos_phonon: [
        {
          energies: [0, 1e-19],
          total: [
            {
              value: [0, 1e18],
              normalization_factor: 1e-19
            }
          ]
        }
      ],
      band_structure_electronic: [
        {
          reciprocal_cell: [[1e9, 0, 0], [0, 1e9, 0], [0, 0, 1e9]],
          segment: [
            {
              energies: [[[0], [1e-19]]],
              kpoints: [[0, 0, 0], [0.5, 0.5, 0.5]],
              endpoints_labels: ['L', 'K']
            }
          ],
          band_gap: [{
            energy_highest_occupied: 0,
            value: 1e-19,
            type: 'indirect'
          }]
        }
      ],
      band_structure_phonon: [
        {
          reciprocal_cell: [[1e9, 0, 0], [0, 1e9, 0], [0, 0, 1e9]],
          segment: [
            {
              energies: [[[0], [1e-19]]],
              kpoints: [[0, 0, 0], [0.5, 0.5, 0.5]],
              endpoints_labels: ['L', 'K']
            }
          ]
        }
      ]
    }
  ]
}
]

// Results for a repository API query
export function getIndex() {
  return {
    ...cloneDeep(common),
    results: cloneDeep({...resultsDftBulk})
  }
}

// Result for an archive API query
export function getArchive() {
  return {
    metadata: {
      ...cloneDeep(common)
    },
    workflow: {...cloneDeep(workflow)},
    results: {...cloneDeep(resultsDftBulk)},
    run: cloneDeep(run)
  }
}
