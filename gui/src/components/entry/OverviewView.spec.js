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
import { waitFor, within } from '@testing-library/dom'
import {render, screen, expectQuantity, readArchive, startAPI, closeAPI} from '../conftest.spec'
import { expectPlotButtons } from '../visualization/conftest.spec'
import {
  expectComposition,
  expectSymmetry,
  expectLatticeParameters
} from './conftest.spec'
import OverviewView from './OverviewView'
import EntryContext from './EntryContext'

test('correctly renders metadata and all properties', async () => {
  await startAPI('tests.states.entry.dft', 'tests/data/entry/dft')
  render(<EntryContext entryId={'dft_bulk'}>
    <OverviewView />
  </EntryContext>)

  // Wait to load the entry metadata, i.e. wait for some of the text to appear
  await screen.findByText('VASP')

  // We read the JSON archive corresponding to the tested API entry. Using this
  // data makes writing assertions much easier.
  const {index} = (await readArchive('../../../tests/states/archives/dft.json'))

  // Check if all method quantities are shown (on the left)
  expectQuantity('results.method.simulation.program_name', 'VASP')
  expectQuantity('results.method.simulation.program_version', '1')
  expectQuantity('results.method.simulation.dft.xc_functional_type', 'GGA')
  expectQuantity('results.method.simulation.dft.xc_functional_names', 'GGA_C_PBE, GGA_X_PBE')
  expectQuantity('results.method.simulation.dft.basis_set_type', 'plane waves')
  expectQuantity('results.method.simulation.dft.basis_set_name', 'STO-3G')
  expectQuantity('results.method.simulation.dft.van_der_Waals_method', 'G06')
  expectQuantity('results.method.simulation.dft.relativity_method', 'scalar_relativistic_atomic_ZORA')

  // Check if all metadata is shown (on the left)
  expectQuantity('results.method.method_name', index)
  expectQuantity('comment', index)
  expectQuantity('references', index.references[0])
  expectQuantity('authors', 'Markus Scheidgen')
  expectQuantity('mainfile', index)
  expectQuantity('entry_id', index)
  expectQuantity('upload_id', index)
  expectQuantity('results.material.material_id', index)
  expectQuantity(undefined, `${index.nomad_version}/${index.nomad_commit}`, 'processing version', 'Version used in the last processing')
  // TODO: add the following to the state for testing.
  // expectQuantity('datasets', index.datasets[0].dataset_name)
  // expectQuantity('upload_create_time', new Date(index.upload_create_time).toLocaleString())
  // expectQuantity('last_processing_time', new Date(index.last_processing_time).toLocaleString())

  // Check if all material data is shown (on the right, in the materials card)
  expectComposition(index)
  expectSymmetry(index)
  expectLatticeParameters(index)
  // expectStructure(index) // TODO: The click introduced here breaks the subsequent tests

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

  // Check if all placeholders disappear
  const dosPhononPlaceholder = screen.queryByTestId('dos-phonon-placeholder')
  const bsPhononPlaceholder = screen.queryByTestId('bs-phonon-placeholder')
  const heatCapacityPlaceholder = screen.queryByTestId('heat-capacity-placeholder')
  const energyFreePlaceholder = screen.queryByTestId('energy-free-placeholder')
  const dosElectronicPlaceholder = screen.queryByTestId('dos-electronic-placeholder')
  const bsElectronicPlaceholder = screen.queryByTestId('bs-electronic-placeholder')
  const energyVolumeCurvePlaceholder = screen.queryByTestId('energy-volume-curve-placeholder')
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

  closeAPI()
})

function expectQuantityToBe(name, label, value, root = screen) {
  const element = root.queryByTestId(`quantity-${name}`)
  expect(within(element).getByText(label)).toBeInTheDocument()
  if (value === undefined) return
  const values = Array.isArray(value) ? value : [value]
  values.forEach(expectedValue => expect(within(element).queryAllByText(expectedValue).length).not.toBe(0))
}

test('eln overview as a reviewer', async () => {
  await startAPI('tests.states.entry.eln', 'tests/data/entry/eln-reviewer', 'ttester', 'password')
  render(<EntryContext entryId={'bC7byHvWJp62Sn9uiuJUB38MT5j-'}>
    <OverviewView />
  </EntryContext>)

  await screen.findByText('HotplateAnnealing')

  expect(screen.queryByTitle("Replace this entry's mainfile")).not.toBeInTheDocument()
  expect(screen.queryByTitle('Save archive')).not.toBeInTheDocument()
  expect(screen.queryByTitle('Delete archive')).not.toBeInTheDocument()

  const sectionCards = screen.queryAllByTestId('property-card')
  expect(sectionCards.length).toBe(3)

  const cardSample = sectionCards[0]
  const cardPvdEvaporation = sectionCards[1]
  const cardHotplateAnnealing = sectionCards[2]

  expect(within(cardSample).getByText('Sample')).toBeVisible()
  expectQuantityToBe('chemical_formula', 'chemical formula', undefined, within(cardSample))
  expectQuantityToBe('name', 'name', 'ELN example sample', within(cardSample))
  expectQuantityToBe('lab_id', 'lab id', '001', within(cardSample))
  expectQuantityToBe('description', 'description', undefined, within(cardSample))
  expectQuantityToBe('tags', 'tags', 'project', within(cardSample))
  expectQuantityToBe('chemicals', 'chemicals', ['../upload/raw/Copper_II_Selenide.archive.json#data', '../upload/raw/Tin_II_Selenide.archive.json#data', '../upload/raw/Zinc_Selenide.archive.json#data'], within(cardSample))
  expectQuantityToBe('substrate_type', 'substrate type', 'SLG', within(cardSample))
  expectQuantityToBe('substrate_thickness', 'substrate thickness', undefined, within(cardSample))
  expectQuantityToBe('sample_is_from_collaboration', 'sample is from collaboration', undefined, within(cardSample))

  expect(within(cardPvdEvaporation).getByText('PvdEvaporation')).toBeVisible()
  expectQuantityToBe('instrument', 'instrument', '../upload/raw/PVD-P.archive.json#data', within(cardPvdEvaporation))
  expectQuantityToBe('row_refs', 'row refs', undefined, within(cardPvdEvaporation))
  expectQuantityToBe('data_file', 'data file', 'PVDProcess.csv', within(cardPvdEvaporation))
  expectQuantityToBe('time', 'time', ['0', '1', '2', '3', '4', 'and 9642 more items'], within(cardPvdEvaporation))
  expectQuantityToBe('chamber_pressure', 'chamber pressure', ['0.00313', '0.00315', '0.00313', '0.00313', '0.00314', 'and 9642 more items'], within(cardPvdEvaporation))
  expectQuantityToBe('substrate_temperature', 'substrate temperature', ['32.4132', '32.4141', '32.416', '32.4175', '32.4181', 'and 9642 more items'], within(cardPvdEvaporation))

  // Test if the plot is there
  expect(within(cardPvdEvaporation).getByText(/Chamber Pressure \(GPa\)/)).toBeVisible()
  expect(within(cardPvdEvaporation).getByText(/Substrate Temperature \(K\)/)).toBeVisible()
  expect(within(cardPvdEvaporation).getByText(/Time \(fs\)/)).toBeVisible()

  expect(within(cardHotplateAnnealing).getByText('HotplateAnnealing')).toBeVisible()
  expectQuantityToBe('instrument', 'instrument', undefined, within(cardHotplateAnnealing))
  expectQuantityToBe('method', 'method', undefined, within(cardHotplateAnnealing))
  expectQuantityToBe('row_refs', 'row refs', undefined, within(cardHotplateAnnealing))
  expectQuantityToBe('set_temperature', 'set temperature', '373.15', within(cardHotplateAnnealing))
  expectQuantityToBe('duration', 'duration', '60', within(cardHotplateAnnealing))

  closeAPI()
})

function expectNumberEditQuantity(numberField, unitField, value, unit) {
  const numberFieldValueInput = within(numberField).getByRole('textbox')
  const numberFieldUnitInput = within(unitField).getByRole('textbox', { hidden: true })
  expect(numberFieldValueInput.value).toEqual(value)
  expect(numberFieldUnitInput.value).toEqual(unit)
}

test.each([
  [
    'an author',
    'tests.states.entry.eln',
    'tests/data/entry/eln-author',
    'bC7byHvWJp62Sn9uiuJUB38MT5j-',
    'test',
    'password'
  ], [
    'a coauthor',
    'tests.states.entry.eln',
    'tests/data/entry/eln-coauthor',
    'bC7byHvWJp62Sn9uiuJUB38MT5j-',
    'scooper',
    'password'
  ]
])('eln overview as %s', async (name, state, snapshot, entryId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<EntryContext entryId={entryId}>
    <OverviewView />
  </EntryContext>)

  await screen.findByText('HotplateAnnealing')

  const saveButton = screen.queryByTitle('Save archive').closest('button')
  expect(saveButton).toBeInTheDocument()
  expect(saveButton).toBeDisabled()

  const reUploadButton = screen.queryByTitle("Replace this entry's mainfile").closest('button')
  expect(reUploadButton).toBeInTheDocument()
  expect(reUploadButton).toBeEnabled()

  const deleteButton = screen.queryByTitle('Delete archive').closest('button')
  expect(deleteButton).toBeInTheDocument()
  expect(deleteButton).toBeEnabled()

  const sectionCards = screen.queryAllByTestId('property-card')
  expect(sectionCards.length).toBe(3)

  const cardSample = sectionCards[0]
  const cardPvdEvaporation = sectionCards[1]
  const cardHotplateAnnealing = sectionCards[2]

  expect(within(cardSample).getByText('Sample')).toBeVisible()
  let numberFieldValue = within(cardSample).queryAllByTestId('number-edit-quantity-value')
  let numberFieldUnit = within(cardSample).queryAllByTestId('number-edit-quantity-unit')
  expectNumberEditQuantity(numberFieldValue[0], numberFieldUnit[0], '', 'Å')

  // Test if the plot is there
  expect(within(cardPvdEvaporation).getByText(/Chamber Pressure \(GPa\)/)).toBeVisible()
  expect(within(cardPvdEvaporation).getByText(/Substrate Temperature \(K\)/)).toBeVisible()
  expect(within(cardPvdEvaporation).getByText(/Time \(fs\)/)).toBeVisible()

  numberFieldValue = within(cardHotplateAnnealing).queryAllByTestId('number-edit-quantity-value')
  numberFieldUnit = within(cardHotplateAnnealing).queryAllByTestId('number-edit-quantity-unit')
  expectNumberEditQuantity(numberFieldValue[0], numberFieldUnit[0], '373.15', 'K')

  closeAPI()
})
