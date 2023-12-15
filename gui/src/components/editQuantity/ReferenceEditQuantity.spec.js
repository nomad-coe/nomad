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
import {render, screen, startAPI, closeAPI, waitForGUI} from '../conftest.spec'
import OverviewView from '../entry/OverviewView'
import { EntryContext } from '../entry/EntryContext'
import userEvent from '@testing-library/user-event'
import {fireEvent} from '@testing-library/react'
import {act} from 'react-dom/test-utils'

const testSectionSelectDialog = async () => {
  const dialog = screen.getByTestId('section-select-dialog')
  expect(within(dialog).queryByText('Filters')).toBeInTheDocument()

  await waitFor(() => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(5))
  const rows = within(dialog).getAllByTestId('datatable-row')

  await waitFor(() => expect(within(rows[0]).queryByText('ref2.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[0]).queryByText('Substance')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[1]).queryByText('ref4.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[1]).queryByText('Structure')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[2]).queryByText('ref3.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[2]).queryByText('Sample')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[3]).queryByText('ref5.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[3]).queryByText('SubstanceList')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[4]).queryByText('SubstanceExtended1 (./data/contents/3)')).toBeInTheDocument())

  const sections = within(rows[4]).queryAllByRole('menuitem')
  expect(sections.length).toBe(4)
  await waitFor(() => expect(within(sections[0]).queryByText('Substance1 (./data/contents/0)')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[0]).queryByTestId('check-icon')).not.toBeInTheDocument())
  await waitFor(() => expect(within(sections[1]).queryByText('AContents1 (./data/contents/1)')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[1]).queryByTestId('check-icon')).not.toBeInTheDocument())
  await waitFor(() => expect(within(sections[2]).queryByText('BContents1 (./data/contents/2)')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[2]).queryByTestId('check-icon')).not.toBeInTheDocument())
  await waitFor(() => expect(within(sections[3]).queryByText('SubstanceExtended1 (./data/contents/3)')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[3]).queryByTestId('check-icon')).toBeInTheDocument())

  // close the dialog
  await userEvent.click(within(dialog).getByRole('button', {name: /cancel/i}))
  await waitFor(() => expect(screen.queryByTestId('section-select-dialog')).not.toBeInTheDocument())
}

const testSectionSelectAutocomplete = async () => {
  await waitFor(() => expect(screen.queryAllByTestId('section-select-entry-deactivate').length).toBe(4))

  const sectionSelectEntries = screen.getAllByTestId('section-select-entry-activated')
  expect(sectionSelectEntries.length).toBe(4)

  await waitFor(() => expect(within(sectionSelectEntries[0]).queryByText('ref2.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[0]).queryByText('upload id: references_upload_id1')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[1]).queryByText('ref3.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[1]).queryByText('upload id: references_upload_id1')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[2]).queryByText('ref4.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[2]).queryByText('upload id: references_upload_id1')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[3]).queryByText('ref5.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[3]).queryByText('upload id: references_upload_id1')).toBeInTheDocument())

  const sectionSelectDeactivateEntries = screen.getAllByTestId('section-select-entry-deactivate')

  await waitFor(() => expect(within(sectionSelectDeactivateEntries[0]).queryByText('correct-reference.data.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectDeactivateEntries[1]).queryByText('lost-reference.data.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectDeactivateEntries[2]).queryByText('ref1.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectDeactivateEntries[3]).queryByText('ref6.archive.yaml')).toBeInTheDocument())

  const sectionSelectPaths = screen.getAllByTestId('section-select-path')
  expect(sectionSelectPaths.length).toBe(4)

  await waitFor(() => expect(within(sectionSelectPaths[0]).queryByText('Substance1 (./data/contents/0)')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectPaths[1]).queryByText('AContents1 (./data/contents/1)')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectPaths[2]).queryByText('BContents1 (./data/contents/2)')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectPaths[3]).queryByText('SubstanceExtended1 (./data/contents/3)')).toBeInTheDocument())

  await act(async () => { userEvent.click(sectionSelectPaths[1]) })
}

const testCreateReferenceDialog = async () => {
  const dialog = screen.getByTestId('create-reference-dialog')
  expect(within(dialog).queryByText('Create new reference of type Substance')).toBeInTheDocument()

  const createButton = within(dialog).getByRole('button', { name: /create/i })
  expect(createButton).toBeDisabled()
  const nameTextField = within(dialog).getByTestId('new-reference-name')
  const nameInput = within(nameTextField).getByRole('textbox')
  fireEvent.change(nameInput, { target: { value: 'newReference' } })

  await waitFor(() => expect(within(dialog).queryByText('File name: newReference.archive.json')).toBeInTheDocument())
  await waitFor(() => expect(createButton).toBeEnabled())

  await userEvent.click(createButton)
  await waitFor(() => expect(screen.queryByTestId('create-reference-dialog')).not.toBeInTheDocument())
}

test.each([
  [
    'correct reference',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity',
    '4WgzB6xcTzWB_Xk9UNUH4HB3IRKJ',
    'test',
    'password'
  ]
])('test referenceEditQuantity %s', async (name, state, snapshot, entryId, username, password) => {
  await startAPI(state, snapshot, username, password)
  await render(
    <EntryContext entryId={entryId}>
      <OverviewView />
    </EntryContext>
  )

  await waitForGUI(2000)

  await waitFor(() => expect(screen.getByText('Powder_mixture')).toBeInTheDocument())

  const sectionCards = screen.queryAllByTestId('property-card')
  const powderMixture = sectionCards[0]

  const referenceEditQuantities = within(powderMixture).getAllByTestId('reference-edit-quantity')
  const referenceEditQuantity = referenceEditQuantities[0]
  const inputTextField = within(referenceEditQuantity).getByRole('textbox')
  await waitFor(() => expect(inputTextField.value).toEqual('SubstanceExtended1'))

  // test section select dialog
  const editReferenceButton = within(referenceEditQuantity).getByTitle('Search for the references').closest('button')
  expect(editReferenceButton).toBeEnabled()

  expect(within(referenceEditQuantity).queryByTitle('Create and assign a new reference')).not.toBeInTheDocument()

  // open edit section select dialog
  await act(async () => { userEvent.click(editReferenceButton) })

  await waitForGUI(1000, true)

  await testSectionSelectDialog()

  // test section select autocomplete
  fireEvent.change(inputTextField, { target: { value: 'ref' } })
  await waitForGUI(2000) // the input changes have been debounced

  await testSectionSelectAutocomplete()

  await waitFor(() => expect(inputTextField.value).toEqual('ref5.archive.yaml#data/contents/1'))

  // test delete reference
  const deleteReferenceButton = within(referenceEditQuantity).getByTitle('Clear').closest('button')
  expect(deleteReferenceButton).toBeEnabled()
  await act(async () => { await userEvent.click(deleteReferenceButton) })

  // when the reference is undefined the create reference button should appear
  const createReferenceButton = within(referenceEditQuantity).getByTitle('Create and assign a new reference').closest('button')
  expect(createReferenceButton).toBeEnabled()
  await act(async () => { await userEvent.click(createReferenceButton) })

  await testCreateReferenceDialog()

  await waitForGUI(1000, true)
  await waitFor(() => expect(inputTextField.value).toEqual('newReference.archive.json'))

  await waitFor(() => expect(within(powderMixture).queryByText('The referenced value does not exist anymore')).not.toBeInTheDocument())
  await waitFor(() => expect(within(powderMixture).queryByText('The provided path does not exist')).not.toBeInTheDocument())

  closeAPI()
})

test.each([
  [
    'lost entry',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity-lost-entry',
    'ScmGivaG2TTQTSYjlJGjIaSF_xTn',
    'test',
    'password',
    '',
    'The referenced value does not exist anymore'
  ],
  [
    'wrong upload id',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity-wrong-upload-id',
    'qZrvjM8MQcd1NX0CYF-2sJqjrgKR',
    'test',
    'password',
    '',
    'The referenced value does not exist anymore'
  ],
  [
    'lost path',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity-lost-path',
    '12GPrF13SLmhgKkvjdCbmNKeMfUv',
    'test',
    'password',
    'ref5.archive.yaml#wrong/path',
    'The provided path does not exist'
  ]
])('test referenceEditQuantity %s', async (name, state, snapshot, entryId, username, password, inputValue, error) => {
  await startAPI(state, snapshot, username, password)
  await render(
    <EntryContext entryId={entryId}>
      <OverviewView />
    </EntryContext>
  )

  await waitForGUI(2000)

  await waitFor(() => expect(screen.getByText('Powder_mixture')).toBeInTheDocument())

  const sectionCards = screen.queryAllByTestId('property-card')
  const powderMixture = sectionCards[0]

  within(powderMixture).getAllByTestId('reference-edit-quantity')
  const referenceEditQuantities = within(powderMixture).getAllByTestId('reference-edit-quantity')
  const referenceEditQuantity = referenceEditQuantities[0]
  const inputTextField = within(referenceEditQuantity).getByRole('textbox')
  await waitFor(() => expect(inputTextField.value).toEqual(inputValue))
  await waitFor(() => expect(within(powderMixture).getByText(error)).toBeInTheDocument())

  closeAPI()
})
