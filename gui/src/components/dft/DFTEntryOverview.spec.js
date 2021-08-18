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
import DFTEntryOverview from './DFTEntryOverview'
import {
  repoDftBulk,
  archiveDftBulk
} from '../../../tests/DFTBulk'
import { useApi } from '../apiV1'

jest.mock('../apiV1')

beforeAll(() => {
  useApi.mockReturnValue({
    results: entry_id => wait(archives.get(entry_id))
  })
})

afterAll(() => jest.unmock('../apiV1'))

test('correctly renders method and material data', async () => {
  const entry = repoDftBulk
  const results = entry.results
  // const archive = archiveDftBulk

  renderWithAPIRouter(
    <DFTEntryOverview
      data={repoDftBulk}
    ></DFTEntryOverview>
  )
  // Wait for the component to update to its final state before performing the
  // tests. This prevents the test from removing the component before the API
  // call is finished (leading to "Can't perform a React state update on an
  // unmounted component...") and also removes the warning about the component
  // doing something unexpected ("An update to ...  inside a test was not
  // wrapped in act(...)")
  const dos_electronic_placeholder = screen.getByTestId('dos-electronic-placeholder')
  await waitFor(() => {
    expect(dos_electronic_placeholder).not.toBeInTheDocument()
  })

  const code_name = screen.getByTitle('The name of the used program.')
  expect(within(code_name).getByText('program name')).toBeInTheDocument()
  expect(within(code_name).getByText(results.method.simulation.program_name)).toBeInTheDocument()
  const code_version = screen.getByTitle('The version of the used program.')
  expect(within(code_version).getByText('program version')).toBeInTheDocument()
  expect(within(code_version).getByText(results.method.simulation.program_version)).toBeInTheDocument()
  // const xc_family = screen.getByTitle('The libXC based xc functional classification used in the simulation.')
  // expect(within(xc_family).getByText('xc functional family')).toBeInTheDocument()
  // expect(within(xc_family).getByText(repo.dft.xc_functional)).toBeInTheDocument()
  // const xc_names = screen.getByTitle('The list of libXC functional names that where used in this entry.')
  // expect(within(xc_names).getByText('xc functional names')).toBeInTheDocument()
  // expect(within(xc_names).getByText(repo.dft.xc_functional_names.join(', '))).toBeInTheDocument()
  // const basis_set_type = screen.getByTitle('The used basis set functions.')
  // expect(within(basis_set_type).getByText('basis set type')).toBeInTheDocument()
  // expect(within(basis_set_type).getByText(repo.dft.basis_set)).toBeInTheDocument()
  // const basis_set_name = screen.getByTitle('Unique string identifying the basis set used for the final wavefunctions calculated with XC_method. It might identify a class of basis sets, often matches one of the strings given in any of basis_set_name.')
  // expect(within(basis_set_name).getByText('basis set name')).toBeInTheDocument()
  // expect(within(basis_set_name).getByText(archive.section_run[0].section_method[0].basis_set)).toBeInTheDocument()
  // const vdw = screen.getByTitle('The used Van der Waals method.')
  // expect(within(vdw).getByText('van der Waals method')).toBeInTheDocument()
  // expect(within(vdw).getByText(archive.section_run[0].section_method[0].van_der_Waals_method)).toBeInTheDocument()
  // const relativity = screen.getByTitle('Describes the relativistic treatment used for the calculation of the final energy and related quantities. If skipped or empty, no relativistic treatment is applied.')
  // expect(within(relativity).getByText('relativity method')).toBeInTheDocument()
  // expect(within(relativity).getByText(archive.section_run[0].section_method[0].relativity_method)).toBeInTheDocument()

  const comment = screen.getByTitle('A user provided comment for this entry')
  expect(within(comment).getByText('comment')).toBeInTheDocument()
  expect(within(comment).getByText(entry.comment)).toBeInTheDocument()
  const references = screen.getByTitle('User provided references (URLs) for this entry')
  expect(within(references).getByText('references')).toBeInTheDocument()
  expect(within(references).getByText(entry.references[0])).toBeInTheDocument()
  const authors = screen.getByTitle('All authors (uploader and co-authors)')
  expect(within(authors).getByText('authors')).toBeInTheDocument()
  expect(within(authors).getByText(entry.authors[0].name)).toBeInTheDocument()
  const datasets = screen.getByTitle('A list of user curated datasets this entry belongs to.')
  expect(within(datasets).getByText('datasets')).toBeInTheDocument()
  expect(within(datasets).getByText(entry.datasets[0].name)).toBeInTheDocument()
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
  const upload_time = screen.getByTitle('The date and time this entry was uploaded to nomad')
  expect(within(upload_time).getByText('upload time')).toBeInTheDocument()
  expect(within(upload_time).getByText(new Date(entry.upload_time).toLocaleString())).toBeInTheDocument()
  const last_processing = screen.getByTitle('The datetime of the last processing')
  expect(within(last_processing).getByText('last processing')).toBeInTheDocument()
  expect(within(last_processing).getByText(new Date(entry.last_processing).toLocaleString())).toBeInTheDocument()
  const processing_version = screen.getByTitle('Version used in the last processing')
  expect(within(processing_version).getByText('processing version')).toBeInTheDocument()
  expect(within(processing_version).getByText(`${entry.nomad_version}/${entry.nomad_commit}`)).toBeInTheDocument()

  const material = results.material
  const formula = screen.getByTitle('The chemical formula for a structure in Hill form with element symbols followed by integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.')
  expect(within(formula).getByText('formula')).toBeInTheDocument()
  expect(within(formula).getByText(material.chemical_formula_hill)).toBeInTheDocument()
  const structural_type = screen.getByTitle('Classification based on structural features.')
  expect(within(structural_type).getByText('structural type')).toBeInTheDocument()
  expect(within(structural_type).getByText(material.structural_type)).toBeInTheDocument()
  const material_name = screen.getByTitle('Meaningful names for this a material if any can be assigned.')
  expect(within(material_name).getByText('material name')).toBeInTheDocument()
  expect(within(material_name).getByText(material.material_name)).toBeInTheDocument()
  const crystal_system = screen.getByTitle('Name of the crystal system. Can be one of the following: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal or cubic.')
  expect(within(crystal_system).getByText('crystal system')).toBeInTheDocument()
  expect(within(crystal_system).getByText(material.symmetry.crystal_system)).toBeInTheDocument()
  const space_group = screen.getByTitle('Space group symbol and number')
  expect(within(space_group).getByText('space group')).toBeInTheDocument()
  expect(within(space_group).getByText(`${material.symmetry.space_group_symbol} (${material.symmetry.space_group_number})`)).toBeInTheDocument()
  const method_name = screen.getByTitle('Common name for the used method.')
  expect(within(method_name).getByText('method name')).toBeInTheDocument()
  expect(within(method_name).getByText('DFT')).toBeInTheDocument()

  // The test DOM does not support canvas or WebGL, and trying to add mocks for
  // them does not seem to work ATM.  Thus we expect a message saying that the
  // 3D viewers are disabled.
  const msgs = screen.getAllByText('Could not display the visualization as your browser does not support WebGL content.')
  expect(msgs).toHaveLength(2)
})

function expectPlotButtons(plot) {
  expect(within(plot).getByRole('button', {name: 'Reset view'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'Toggle fullscreen'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'Capture image'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'View data in the archive'})).toBeInTheDocument()
}

test('correctly renders electronic properties', async () => {
  renderWithAPIRouter(
    <DFTEntryOverview
      data={repoDftBulk}
    ></DFTEntryOverview>
  )

  expect(screen.getByText('Electronic properties')).toBeInTheDocument()
  expect(screen.getByText('Band structure')).toBeInTheDocument()
  expect(screen.getByText('Density of states')).toBeInTheDocument()
  expect(screen.getByText('Brillouin zone')).toBeInTheDocument()

  // Placeholders are shown in the beginning but removed when plot is loaded.
  const dosElectronicPlaceholder = screen.getByTestId('dos-electronic-placeholder')
  expect(dosElectronicPlaceholder).toBeInTheDocument()
  const bsElectronicPlaceholder = screen.getByTestId('bs-electronic-placeholder')
  expect(bsElectronicPlaceholder).toBeInTheDocument()
  await waitFor(() => { expect(dosElectronicPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(bsElectronicPlaceholder).not.toBeInTheDocument() })

  // Check that plot buttons are displayed
  const dosElectronic = screen.getByTestId('dos-electronic')
  expectPlotButtons(dosElectronic)
  const bsElectronic = screen.getByTestId('bs-electronic')
  expectPlotButtons(bsElectronic)
})

test('correctly renders vibrational properties', async () => {
  renderWithAPIRouter(
    <DFTEntryOverview
      data={repoDftBulk}
    ></DFTEntryOverview>
  )

  expect(screen.getByText('Vibrational properties')).toBeInTheDocument()
  expect(screen.getByText('Phonon dispersion')).toBeInTheDocument()
  expect(screen.getByText('Phonon density of states')).toBeInTheDocument()
  expect(screen.getByText('Heat capacity')).toBeInTheDocument()
  expect(screen.getByText('Helmholtz free energy')).toBeInTheDocument()

  // Placeholders are shown in the beginning but removed when plot is loaded.
  const dosPhononPlaceholder = screen.getByTestId('dos-phonon-placeholder')
  expect(dosPhononPlaceholder).toBeInTheDocument()
  const bsPhononPlaceholder = screen.getByTestId('bs-phonon-placeholder')
  expect(bsPhononPlaceholder).toBeInTheDocument()
  const heatCapacityPlaceholder = screen.getByTestId('heat-capacity-placeholder')
  expect(heatCapacityPlaceholder).toBeInTheDocument()
  const energyFreePlaceholder = screen.getByTestId('energy-free-placeholder')
  expect(energyFreePlaceholder).toBeInTheDocument()
  await waitFor(() => { expect(dosPhononPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(bsPhononPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(heatCapacityPlaceholder).not.toBeInTheDocument() })
  await waitFor(() => { expect(energyFreePlaceholder).not.toBeInTheDocument() })

  // Check that plot buttons are displayed
  const bsPhonon = screen.getByTestId('bs-phonon')
  expectPlotButtons(bsPhonon)
  const dosPhonon = screen.getByTestId('dos-phonon')
  expectPlotButtons(dosPhonon)
  const heatCapacity = screen.getByTestId('heat-capacity')
  expectPlotButtons(heatCapacity)
  const energyFree = screen.getByTestId('energy-free')
  expectPlotButtons(energyFree)
})
