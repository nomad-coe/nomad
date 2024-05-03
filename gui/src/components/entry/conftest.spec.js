
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

import { within } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { expectQuantity, screen } from '../conftest.spec'
import { expectPlotButtons } from '../visualization/conftest.spec'
import { traverseDeep, serializeMetainfo } from '../../utils'
import { ui } from '../../config'

/*****************************************************************************/
// Expects
export function expectComposition(index) {
  expectQuantity(undefined, index.results.material.chemical_formula_hill, 'Formula', 'The chemical formula that describes the simulated system or experimental sample.')
  expectQuantity('results.material.structural_type', index)
  expectQuantity('results.material.elements', index)
  expectQuantity('results.material.n_elements', '1 (unary)', 'N elements')
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

export async function expectStructure(index, root = screen) {
  // Before clicking only the default value should be visible
  const select = root.getByText('Original')
  expect(select).toBeInTheDocument()
  const conv = root.queryByText('Conventional')
  expect(conv).not.toBeInTheDocument()
  const prim = root.queryByText('Primitive')
  expect(prim).not.toBeInTheDocument()

  // After clicking all the options should be shown
  await userEvent.click(select)
  expect(await root.findByText('Conventional')).toBeInTheDocument()
  expect(await root.findByText('Primitive')).toBeInTheDocument()
  await userEvent.click(select)
}

export function expectNoStructure(index, root = screen) {
  const select = root.getByText('Structure')
  const parent = select.parentElement
  expect(within(parent).getByText('no data')).toBeInTheDocument()
}

export function expectTrajectory(index, root = screen) {
  expect(root.queryByText('Trajectory')).toBeInTheDocument()
  expectQuantity('results.properties.thermodynamic.trajectory.provenance.molecular_dynamics.time_step', '1 fs')
  expectQuantity('results.properties.thermodynamic.trajectory.provenance.molecular_dynamics.ensemble_type', 'NVT')
  const trajectoryPlot = screen.getByTestId('trajectory')
  expectPlotButtons(trajectoryPlot)
}

/**
 * Tests that an methodology itme is displayed correctly.
 *
 * @param {str} title The methodology title.
 * @param {object} data The data that should be shown.
 * @param {str} path The path for the methodology information.
 * @param {object} container The root element to perform the search on.
 */
export async function expectMethodologyItem(
  title,
  data,
  path,
  container = document.body
) {
  const root = within(container)
  if (!data) {
      expect(root.queryByText(title)).not.toBeInTheDocument()
  } else {
    expect(root.getByText(title)).toBeInTheDocument()
    for (const [key, value] of traverseDeep(data, true)) {
      const quantity = `${path}.${key.join('.')}`
      expectQuantity(quantity, serializeMetainfo(quantity, value, ui.unit_systems.options.Custom.units))
    }
  }
}
