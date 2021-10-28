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
import 'regenerator-runtime/runtime'
import { renderWithAPIRouter, archives, wait } from '../../testutils'
import { screen } from '@testing-library/react'
import { waitFor, within } from '@testing-library/dom'
import '@testing-library/jest-dom/extend-expect'
import {
  repoDftBulk
} from '../../../tests/DFTBulk'
import { useApi } from '../api'
import OverviewView from './OverviewView'

jest.mock('../api')
const { ResizeObserver } = window

beforeAll(() => {
  // API mock init
  useApi.mockReturnValue({
    api: {
      results: entry_id => wait(archives.get(entry_id)),
      entry: () => wait({
        entry_id: repoDftBulk.entry_id,
        data: repoDftBulk
      })
    }
  })

  // ResizeObserver mock init
  delete window.ResizeObserver
  window.ResizeObserver = jest.fn().mockImplementation(() => ({
    observe: jest.fn(),
    unobserve: jest.fn(),
    disconnect: jest.fn()
  }))
})

afterAll(() => {
  // API mock cleanup
  jest.unmock('../api')

  // ResizeObserver mock cleanup
  window.ResizeObserver = ResizeObserver
  jest.restoreAllMocks()
})

function expectPlotButtons(plot) {
  expect(within(plot).getByRole('button', {name: 'Reset view'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'Toggle fullscreen'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'Capture image'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'View data in the archive'})).toBeInTheDocument()
}

test('correctly renders metadata and all properties', async () => {
  const entry = repoDftBulk
  const results = entry.results

  renderWithAPIRouter(
    <OverviewView entryId={repoDftBulk.entry_id} uploadId="dont-care" />
  )

  // Wait to load the entry metadata, i.e. wait for an arbitrary placeholder to appear
  await waitFor(() => {
    expect(screen.getByTestId('dos-electronic-placeholder')).toBeInTheDocument()
  })

  // Check if all method quantities are shown (on the left)
  const method = results.method.simulation
  const code_name = screen.getByTitle('The name of the used program.')
  expect(within(code_name).getByText('program name')).toBeInTheDocument()
  expect(within(code_name).getByText(method.program_name)).toBeInTheDocument()
  const code_version = screen.getByTitle('The version of the used program.')
  expect(within(code_version).getByText('program version')).toBeInTheDocument()
  expect(within(code_version).getByText(method.program_version)).toBeInTheDocument()
  const xc_family = screen.getByTitle('The libXC based xc functional classification used in the simulation.')
  expect(within(xc_family).getByText('xc functional type')).toBeInTheDocument()
  expect(within(xc_family).getByText(method.dft.xc_functional_type)).toBeInTheDocument()
  const xc_names = screen.getByTitle('The list of libXC functional names that where used in this entry.')
  expect(within(xc_names).getByText('xc functional names')).toBeInTheDocument()
  expect(within(xc_names).getByText(method.dft.xc_functional_names.join(', '))).toBeInTheDocument()
  const basis_set_type = screen.getByTitle('The used basis set functions.')
  expect(within(basis_set_type).getByText('basis set type')).toBeInTheDocument()
  expect(within(basis_set_type).getByText(method.dft.basis_set_type)).toBeInTheDocument()
  const basis_set_name = screen.getByTitle('Identifies the basis set.')
  expect(within(basis_set_name).getByText('basis set name')).toBeInTheDocument()
  expect(within(basis_set_name).getByText(method.dft.basis_set_name)).toBeInTheDocument()
  const vdw = screen.getByTitle('The used van der Waals method.')
  expect(within(vdw).getByText('van der Waals method')).toBeInTheDocument()
  expect(within(vdw).getByText(method.dft.van_der_Waals_method)).toBeInTheDocument()
  const relativity = screen.getByTitle('Describes the relativistic treatment used for the calculation of the final energy and related quantities. If skipped or empty, no relativistic treatment is applied.')
  expect(within(relativity).getByText('relativity method')).toBeInTheDocument()
  expect(within(relativity).getByText(method.dft.relativity_method)).toBeInTheDocument()

  // Check if all metadata is shown (on the left)
  const comment = screen.getByTitle('A user provided comment for this entry')
  expect(within(comment).getByText('comment')).toBeInTheDocument()
  expect(within(comment).getByText(entry.comment)).toBeInTheDocument()
  const references = screen.getByTitle('User provided references (URLs) for this entry')
  expect(within(references).getByText('references')).toBeInTheDocument()
  expect(within(references).getByText(entry.references[0])).toBeInTheDocument()
  const authors = screen.getByTitle('All authors (main author and co-authors)')
  expect(within(authors).getByText('authors')).toBeInTheDocument()
  expect(within(authors).getByText(entry.authors[0].name)).toBeInTheDocument()
  const datasets = screen.getByTitle('A list of user curated datasets this entry belongs to.')
  expect(within(datasets).getByText('datasets')).toBeInTheDocument()
  expect(within(datasets).getByText(entry.datasets[0].dataset_name)).toBeInTheDocument()
  const mainfile = screen.getByTitle('The path to the mainfile from the root directory of the uploaded files')
  expect(within(mainfile).getByText('mainfile')).toBeInTheDocument()
  expect(within(mainfile).getByText(entry.mainfile)).toBeInTheDocument()
  const entry_id = screen.getByTitle('The unique primary id for this entry.')
  expect(within(entry_id).getByText('entry id')).toBeInTheDocument()
  expect(within(entry_id).getByText(entry.entry_id)).toBeInTheDocument()
  const material_id = screen.getByTitle('A fixed length, unique material identifier in the form of a hash digest.')
  expect(within(material_id).getByText('material id')).toBeInTheDocument()
  expect(within(material_id).getByText(results.material.material_id)).toBeInTheDocument()
  const upload_id = screen.getByTitle('The persistent and globally unique identifier for the upload of the entry')
  expect(within(upload_id).getByText('upload id')).toBeInTheDocument()
  expect(within(upload_id).getByText(entry.upload_id)).toBeInTheDocument()
  const upload_create_time = screen.getByTitle('The date and time when the upload was created in nomad')
  expect(within(upload_create_time).getByText('upload time')).toBeInTheDocument()
  expect(within(upload_create_time).getByText(new Date(entry.upload_create_time).toLocaleString())).toBeInTheDocument()
  const last_processing_time = screen.getByTitle('The date and time of the last processing.')
  expect(within(last_processing_time).getByText('last processing time')).toBeInTheDocument()
  expect(within(last_processing_time).getByText(new Date(entry.last_processing_time).toLocaleString())).toBeInTheDocument()
  const processing_version = screen.getByTitle('Version used in the last processing')
  expect(within(processing_version).getByText('processing version')).toBeInTheDocument()
  expect(within(processing_version).getByText(`${entry.nomad_version}/${entry.nomad_commit}`)).toBeInTheDocument()

  // Check if all material data is shown (on the right, in the materials card)
  const material = results.material
  const formula = screen.getByTitle('The chemical formula that describes the simulated system or experiment sample.')
  expect(within(formula).getByText('formula')).toBeInTheDocument()
  expect(within(formula).getByText(material.chemical_formula_hill)).toBeInTheDocument()
  const structural_type = screen.getByTitle('Classification based on structural features.')
  expect(within(structural_type).getByText('structural type')).toBeInTheDocument()
  expect(within(structural_type).getByText(material.structural_type)).toBeInTheDocument()
  const material_name = screen.getByTitle('Meaningful names for this a material if any can be assigned.')
  expect(within(material_name).getByText('material name')).toBeInTheDocument()
  expect(within(material_name).getByText(material.material_name)).toBeInTheDocument()
  const crystal_system = screen.getByTitle('Name of the crystal system.')
  expect(within(crystal_system).getByText('crystal system')).toBeInTheDocument()
  expect(within(crystal_system).getByText(material.symmetry.crystal_system)).toBeInTheDocument()
  const space_group = screen.getByTitle('Space group symbol and number')
  expect(within(space_group).getByText('space group')).toBeInTheDocument()
  expect(within(space_group).getByText(`${material.symmetry.space_group_symbol} (${material.symmetry.space_group_number})`)).toBeInTheDocument()
  const method_name = screen.getByTitle('Common name for the used method.')
  expect(within(method_name).getByText('method name')).toBeInTheDocument()
  expect(within(method_name).getByText('DFT')).toBeInTheDocument()

  // Check if all the property cards are shown
  expect(screen.getByText('Electronic properties')).toBeInTheDocument()
  expect(screen.getByText('Band structure')).toBeInTheDocument()
  expect(screen.getByText('Density of states')).toBeInTheDocument()
  expect(screen.getByText('Brillouin zone')).toBeInTheDocument()

  expect(screen.getByText('Vibrational properties')).toBeInTheDocument()
  expect(screen.getByText('Phonon dispersion')).toBeInTheDocument()
  expect(screen.getByText('Phonon density of states')).toBeInTheDocument()
  expect(screen.getByText('Heat capacity')).toBeInTheDocument()
  expect(screen.getByText('Helmholtz free energy')).toBeInTheDocument()
  expect(screen.getByText('Mechanical properties')).toBeInTheDocument()
  expect(screen.getByText('Energy-volume curve')).toBeInTheDocument()
  expect(screen.getByText('Bulk modulus')).toBeInTheDocument()
  expect(screen.getByText('Shear modulus')).toBeInTheDocument()

  // Check if all visualization placeholders are shown
  const dosPhononPlaceholder = screen.getByTestId('dos-phonon-placeholder')
  expect(dosPhononPlaceholder).toBeInTheDocument()
  const bsPhononPlaceholder = screen.getByTestId('bs-phonon-placeholder')
  expect(bsPhononPlaceholder).toBeInTheDocument()
  const heatCapacityPlaceholder = screen.getByTestId('heat-capacity-placeholder')
  expect(heatCapacityPlaceholder).toBeInTheDocument()
  const energyFreePlaceholder = screen.getByTestId('energy-free-placeholder')
  expect(energyFreePlaceholder).toBeInTheDocument()
  const dosElectronicPlaceholder = screen.getByTestId('dos-electronic-placeholder')
  expect(dosElectronicPlaceholder).toBeInTheDocument()
  const bsElectronicPlaceholder = screen.getByTestId('bs-electronic-placeholder')
  expect(bsElectronicPlaceholder).toBeInTheDocument()
  const energyVolumeCurvePlaceholder = screen.getByTestId('energy-volume-curve-placeholder')
  expect(energyVolumeCurvePlaceholder).toBeInTheDocument()

  // Check if all placeholders disappear
  await waitFor(() => { expect(dosElectronicPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(bsElectronicPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(dosPhononPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(bsPhononPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(heatCapacityPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(energyFreePlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(energyVolumeCurvePlaceholder).not.toBeInTheDocument() })

  // The test DOM does not support canvas or WebGL, and trying to add mocks for
  // them does not seem to work ATM. Thus we expect a message saying that the
  // 3D viewers are disabled.
  const msgs = screen.getAllByText('Could not display the visualization as your browser does not support WebGL content.')
  expect(msgs).toHaveLength(2)

  // Check that plot buttons are displayed
  const dosElectronic = screen.getByTestId('dos-electronic')
  expectPlotButtons(dosElectronic)
  const bsElectronic = screen.getByTestId('bs-electronic')
  expectPlotButtons(bsElectronic)
  const bsPhonon = screen.getByTestId('bs-phonon')
  expectPlotButtons(bsPhonon)
  const dosPhonon = screen.getByTestId('dos-phonon')
  expectPlotButtons(dosPhonon)
  const heatCapacity = screen.getByTestId('heat-capacity')
  expectPlotButtons(heatCapacity)
  const energyFree = screen.getByTestId('energy-free')
  expectPlotButtons(energyFree)
  const energyVolumeCurve = screen.getByTestId('energy-volume-curve')
  expectPlotButtons(energyVolumeCurve)

  // Check that tables are shown
  const bulkModulus = screen.getByTestId('bulk-modulus')
  expect(within(bulkModulus).getByText('Type')).toBeInTheDocument()
  expect(within(bulkModulus).getByText('Value (GPa)')).toBeInTheDocument()
  expect(within(bulkModulus).getByText('murnaghan')).toBeInTheDocument()
  const shearModulus = screen.getByTestId('shear-modulus')
  expect(within(shearModulus).getByText('Type')).toBeInTheDocument()
  expect(within(shearModulus).getByText('Value (GPa)')).toBeInTheDocument()
  expect(within(shearModulus).getByText('voigt_reuss_hill_average')).toBeInTheDocument()
})
