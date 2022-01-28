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
import React from 'react'
import { render, screen, within } from '../../../testSetup'
import { expectQuantity } from '../../../testHelpers'
import userEvent from '@testing-library/user-event'
import MaterialCard from './MaterialCard'
import { getIndex, getArchive } from '../../../../tests/DFTBulk'

const index = getIndex()
const archive = getArchive()
const props = new Set(index?.results
  ? index.results.properties.available_properties
  : [])

export function testComposition(index) {
  expectQuantity(undefined, index.results.material.chemical_formula_hill, 'formula', 'The chemical formula that describes the simulated system or experimental sample.')
  expectQuantity('results.material.structural_type', index)
  expectQuantity('results.material.material_name', index)
  expectQuantity('results.material.elements', index)
  expectQuantity('results.material.n_elements', '1 (unary)', 'number of elements')
}

export function testSymmetry(index) {
  expectQuantity('results.material.symmetry.crystal_system', index)
  expectQuantity('results.material.symmetry.bravais_lattice', index)
  expectQuantity('results.material.symmetry.space_group_number', index)
  expectQuantity('results.material.symmetry.space_group_symbol', index)
  expectQuantity('results.material.symmetry.point_group', index)
  expectQuantity('results.material.symmetry.structure_name', index)
}

export function testNoSymmetry(index) {
  const select = screen.getByText('Symmetry')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}

export function testLatticeParameters(index) {
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.a', '5 Å')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.b', '5 Å')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.c', '5 Å')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.alpha', '90 °', 'α')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.beta', '90 °', 'β')
  expectQuantity('results.properties.structures.structure_original.lattice_parameters.gamma', '90 °', 'γ')
  expectQuantity('results.properties.structures.structure_original.cell_volume', '125 Å ^ 3')
}

export function testNoLatticeParameters(index) {
  const select = screen.getByText('Lattice parameters')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}

export function testStructure(index) {
  // Before clicking only the default value should be visible
  const select = screen.getByText('Original')
  expect(select).toBeInTheDocument()
  const conv = screen.queryByText('Conventional')
  expect(conv).not.toBeInTheDocument()
  const prim = screen.queryByText('Primitive')
  expect(prim).not.toBeInTheDocument()

  // After clicking all the options should be shown
  userEvent.click(select)
  expect(screen.queryByText('Conventional')).toBeInTheDocument()
  expect(screen.queryByText('Primitive')).toBeInTheDocument()
  userEvent.click(select)
}

export function testNoStructure(index) {
  const select = screen.getByText('Structure')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}

test('correctly renders entry with all material information', async () => {
  render(<MaterialCard index={index} properties={props} archive={archive}/>)
  testComposition(index)
  testSymmetry(index)
  testLatticeParameters(index)
  testStructure(index)
})

test('correctly renders material without symmetry information', async () => {
  let index = getIndex()
  delete index.results.material.symmetry
  render(<MaterialCard index={index} properties={props} archive={archive}/>)
  testComposition(index)
  testNoSymmetry(index)
  testLatticeParameters(index)
  testStructure(index)
})

test('correctly renders material without lattice information', async () => {
  const index = getIndex()
  delete index.results.properties.structures.structure_original.lattice_parameters
  delete index.results.properties.structures.structure_conventional.lattice_parameters
  delete index.results.properties.structures.structure_primitive.lattice_parameters
  render(<MaterialCard index={index} properties={props} archive={archive}/>)
  testComposition(index)
  testSymmetry(index)
  testNoLatticeParameters(index)
  testStructure(index)
})

test('correctly renders material without any structure information', async () => {
  const index = getIndex()
  delete index.results.properties.structures
  render(<MaterialCard index={index} properties={props} archive={archive}/>)
  testComposition(index)
  testSymmetry(index)
  testNoLatticeParameters(index)
  testNoStructure(index)
})
