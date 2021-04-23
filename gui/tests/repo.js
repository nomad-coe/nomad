/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the 'License');
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an 'AS IS' BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Mock result for the API result of a repository call.
export const repoDftBulk = {
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
  formula: 'Si2',
  dft: {
    code_name: 'VASP',
    code_version: '1',
    xc_functional: 'GGA',
    xc_functional_names: ['GGA_C_PBE', 'GGA_X_PBE'],
    basis_set: 'plane waves',
    system: 'bulk',
    crystal_system: 'cubic',
    spacegroup_symbol: 'Fd-3m',
    spacegroup: 227,
    searchable_quantities: ['electronic_dos']
  },
  encyclopedia: {
    material: {
      material_id: 'Mock material id',
      material_type: 'bulk',
      material_name: 'Silicon'
    }
  },
  results:  {
    material: {
      material_id: 'Mock material id',
      material_name: 'Silicon',
      type_structural: 'bulk',
      chemical_formula_reduced: 'Si2',
      chemical_formula_hill: 'Si2',
      chemical_formula_anonymous: 'A2',
      chemical_formula_descriptive: 'Si2',
      symmetry: {
        crystal_system: 'cubic',
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
          xc_functional_type: 'GGA',
          xc_functional_names: ['GGA_C_PBE', 'GGA_X_PBE']
        }
      }
    },
    properties: {
      available_properties: ["electronic.dos_electronic"],
      electronic: {
        dos_electronic: {
          energy_highest_occupied: [0],
          energy_lowest_unoccupied: [1],
          energies: [0, 1],
          densities: [0, 1],
        }
      }
    }
  }
}
