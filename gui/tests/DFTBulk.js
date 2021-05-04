/*
 * Copyright The NOMAD Authors.
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

// Indexed data that is shared between section_metadata and section_results
const common = {
  domain: 'dft',
  upload_id: 'mock',
  upload_time: '2021-03-17T13:47:32.899000',
  last_processing: '2021-03-17T15:47:32.899000',
  nomad_version: '0.10.0',
  nomad_commit: 'bf3c06fa',
  calc_id: 'dft_bulk',
  entry_id: 'dft_bulk',
  comment: 'Mocked',
  references: ['doi'],
  authors: [{name: 'Lauri Himanen'}],
  datasets: [{dataset_id: 'Mock dataset', name: 'Mock dataset'}],
  mainfile: 'vasp.xml',
  formula: 'Si2'
}
const commonOld = {
  ...common, 
  calc_id: 'dft_bulk_old',
  entry_id: 'dft_bulk_old',
}

const workflow = {
  workflow_type: 'phonon',
  calculation_result_ref: '/section_run/0/section_single_configuration_calculation/0'
}

const materialName = 'Silicon'
const materialType = 'bulk'
const crystalSystem = 'cubic'
const vdwMethod = 'G06'
const relativityMethod = 'scalar_relativistic_atomic_ZORA'
const basisSetName = 'STO-3G'

// Indexed data that is specific to section_metadata
const metaDftBulk = {
  dft: {
    code_name: 'VASP',
    code_version: '1',
    xc_functional: 'GGA',
    xc_functional_names: ['GGA_C_PBE', 'GGA_X_PBE'],
    basis_set: 'plane waves',
    system: materialType,
    crystal_system: crystalSystem,
    spacegroup_symbol: 'Fd-3m',
    spacegroup: 227,
    searchable_quantities: [
      'electronic_dos',
      'electronic_band_structure',
      'phonon_dos',
      'phonon_band_structure',
      'thermodynamical_property_heat_capacity_C_v',
      'vibrational_free_energy_at_constant_volume'
    ]
  },
  encyclopedia: {
    material: {
      material_id: 'Mock material id',
      material_type: materialType,
      material_name: materialName,
      bulk: {
        crystal_system: crystalSystem
      }
    }
  }
}

// Indexed data that is specific to results
const resultsDftBulk = {
  material: {
    material_id: 'Mock material id',
    material_name: materialName,
    type_structural: materialType,
    chemical_formula_reduced: 'Si2',
    chemical_formula_hill: 'Si2',
    chemical_formula_anonymous: 'A2',
    chemical_formula_descriptive: 'Si2',
    symmetry: {
      crystal_system: crystalSystem,
      space_group_symbol: 'Fd-3m',
      space_group_number: 227,
    }
  },
  method: {
    method_name: 'DFT',
    simulation: {
      program_name: 'VASP',
      program_version: '1',
      dft: {
        basis_set_type: 'plane waves',
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
      'energy_free_helmholtz'
    ],
    electronic: {
      dos_electronic: {
        energies: "/section_run/0/section_single_configuration_calculation/0/section_dos/0/dos_energies_normalized",
        densities: "/section_run/0/section_single_configuration_calculation/0/section_dos/0/dos_values_normalized",
        energy_references: "/section_run/0/section_single_configuration_calculation/0/section_dos/0/energy_references"
      },
      band_structure_electronic: {
        segments: ["/section_run/0/section_single_configuration_calculation/0/section_k_band/0/section_k_band_segment/0"],
        reciprocal_cell: "/section_run/0/section_single_configuration_calculation/0/section_k_band/0/reciprocal_cell",
        energy_references: "/section_run/0/section_single_configuration_calculation/0/section_k_band/0/energy_references"
      }
    },
    vibrational: {
      dos_phonon: {
        energies: "/section_run/0/section_single_configuration_calculation/0/section_dos/1/dos_energies",
        densities: "/section_run/0/section_single_configuration_calculation/0/section_dos/1/dos_values_normalized",
      },
      band_structure_phonon: {
        segments: ["/section_run/0/section_single_configuration_calculation/0/section_k_band/1/section_k_band_segment/0"],
      },
      heat_capacity_constant_volume: {
        heat_capacities: "/section_run/0/section_frame_sequence/0/section_thermodynamical_properties/0/thermodynamical_property_heat_capacity_C_v",
        temperatures: "/section_run/0/section_frame_sequence/0/section_thermodynamical_properties/0/thermodynamical_property_temperature"
      },
      energy_free_helmholtz: {
        energies: "/section_run/0/section_frame_sequence/0/section_thermodynamical_properties/0/vibrational_free_energy_at_constant_volume",
        temperatures: "/section_run/0/section_frame_sequence/0/section_thermodynamical_properties/0/thermodynamical_property_temperature"
      }
    }
  }
}

// Section run
const run = [{
    program_name: "VASP",
    section_method: [
      {
        electronic_structure_method: "DFT",
        basis_set: basisSetName,
        van_der_Waals_method: vdwMethod,
        relativity_method: relativityMethod
      }
    ],
    section_frame_sequence: [{
      section_thermodynamical_properties: [{
        thermodynamical_property_heat_capacity_C_v: [0, 1],
        vibrational_free_energy_at_constant_volume: [0, 1],
        thermodynamical_property_temperature: [0, 100]
      }]
    }],
    section_single_configuration_calculation: [
      {
        section_dos: [
          {
            dos_kind: 'electronic',
            dos_energies: [0, 1e-19],
            dos_energies_normalized: [0, 1e-19],
            dos_values: [[0, 1e18]],
            dos_values_normalized: [[0, 1e18]],
            energy_references: [{
              energy_highest_occupied: 0
            }]
          },
          {
            dos_kind: 'vibrational',
            dos_energies: [0, 1e-19],
            dos_energies_normalized: [0, 1e-19],
            dos_values: [[0, 1e18]],
            dos_values_normalized: [[0, 1e18]],
          }
        ],
        section_k_band: [
          {
            band_structure_kind: 'electronic',
            reciprocal_cell: [[1e9, 0, 0], [0, 1e9, 0], [0, 0, 1e9]],
            section_k_band_segment: [
              {
                band_energies: [[[0], [1e-19]]],
                band_k_points: [[0, 0, 0], [0.5, 0.5, 0.5]],
                band_segm_labels: ["L", "K"]
              }
            ],
            energy_references: [{
              energy_highest_occupied: 0
            }]
          },
          {
            band_structure_kind: 'vibrational',
            section_k_band_segment: [
              {
                band_energies: [[[0], [1e-19]]],
                band_k_points: [[0, 0, 0], [0.5, 0.5, 0.5]],
                band_segm_labels: ["L", "K"]
              }
            ]
          }
        ]
      }
    ]
  }
]

// Results for a repository API query
export const repoDftBulk = {
  ...common,
  ...metaDftBulk,
  results: {...resultsDftBulk}
}

export const repoDftBulkOld = {
  ...commonOld,
  ...metaDftBulk
}

// Result for an archive API query
export const archiveDftBulk = {
  section_metadata: {
    ...common,
    ...metaDftBulk,
  },
  section_workflow: {...workflow},
  results: {...resultsDftBulk},
  section_run: run
}

export const archiveDftBulkOld = {
  section_metadata: {
    ...commonOld,
    ...metaDftBulk,
  },
  section_workflow: {...workflow},
  section_run: run
}
