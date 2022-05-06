
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

import { expectQuantity, screen } from '../conftest.spec'
import { within } from '@testing-library/react'
import userEvent from '@testing-library/user-event'

/*****************************************************************************/
// Expects
export function expectComposition(index) {
  expectQuantity(undefined, index.results.material.chemical_formula_hill, 'formula', 'The chemical formula that describes the simulated system or experimental sample.')
  expectQuantity('results.material.structural_type', index)
  expectQuantity('results.material.material_name', index)
  expectQuantity('results.material.elements', index)
  expectQuantity('results.material.n_elements', '1 (unary)', 'number of elements')
}

export function expectSymmetry(index) {
  expectQuantity('results.material.symmetry.crystal_system', index)
  expectQuantity('results.material.symmetry.bravais_lattice', index)
  expectQuantity('results.material.symmetry.space_group_number', index)
  expectQuantity('results.material.symmetry.space_group_symbol', index)
  expectQuantity('results.material.symmetry.point_group', index)
  expectQuantity('results.material.symmetry.structure_name', index)
}

export function expectNoSymmetry(index, root = screen) {
  const select = root.getByText('Symmetry')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}

export function expectLatticeParameters(index) {
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.a', '5 Å')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.b', '5 Å')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.c', '5 Å')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.alpha', '90 °', 'α')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.beta', '90 °', 'β')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.gamma', '90 °', 'γ')
  expectQuantity('results.properties.structures.structure_original.cell_volume', '125 Å^3')
}

export function expectNoLatticeParameters(index, root = screen) {
  const select = root.getByText('Lattice parameters')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}

export function expectStructure(index, root = screen) {
  // Before clicking only the default value should be visible
  const select = root.getByText('Original')
  expect(select).toBeInTheDocument()
  const conv = root.queryByText('Conventional')
  expect(conv).not.toBeInTheDocument()
  const prim = root.queryByText('Primitive')
  expect(prim).not.toBeInTheDocument()

  // After clicking all the options should be shown
  userEvent.click(select)
  expect(root.queryByText('Conventional')).toBeInTheDocument()
  expect(root.queryByText('Primitive')).toBeInTheDocument()
  userEvent.click(select)
}

export function expectNoStructure(index, root = screen) {
  const select = root.getByText('Structure')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}
