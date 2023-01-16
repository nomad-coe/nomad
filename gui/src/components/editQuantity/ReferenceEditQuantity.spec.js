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
import EntryPageContext from '../entry/EntryPageContext'
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

  await waitFor(() => expect(within(rows[2]).queryByText('mySubstance2 (./data/mySample/mySubstance)')).toBeInTheDocument())

  const sections = within(rows[2]).queryAllByRole('menuitem')
  expect(sections.length).toBe(2)
  await waitFor(() => expect(within(sections[0]).queryByText('mySubstance1 (./data/mySubstance)')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[0]).queryByTestId('check-icon')).not.toBeInTheDocument())
  await waitFor(() => expect(within(sections[1]).queryByText('mySubstance2 (./data/mySample/mySubstance)')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[1]).queryByTestId('check-icon')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[3]).queryByText('ref3.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[3]).queryByText('Sample')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[4]).queryByText('ref5.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[4]).queryByText('SubstanceList')).toBeInTheDocument())

  // close the dialog
  await userEvent.click(within(dialog).getByRole('button', {name: /cancel/i}))
  await waitFor(() => expect(screen.queryByTestId('section-select-dialog')).not.toBeInTheDocument())
}

const testSectionSelectAutocomplete = async () => {
  await waitFor(() => expect(screen.queryAllByTestId('section-select-entry-deactivate').length).toBe(4))

  const sectionSelectEntries = screen.getAllByTestId('section-select-entry-activated')
  expect(sectionSelectEntries.length).toBe(4)

  await waitFor(() => expect(within(sectionSelectEntries[0]).queryByText('ref5.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[0]).queryByText('upload id: references_upload_id')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[1]).queryByText('ref3.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[1]).queryByText('upload id: references_upload_id')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[2]).queryByText('ref4.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[2]).queryByText('upload id: references_upload_id')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[3]).queryByText('ref2.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[3]).queryByText('upload id: references_upload_id')).toBeInTheDocument())

  const sectionSelectDeactivateEntries = screen.getAllByTestId('section-select-entry-deactivate')

  await waitFor(() => expect(within(sectionSelectDeactivateEntries[0]).queryByText('ref1.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectDeactivateEntries[1]).queryByText('ref6.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectDeactivateEntries[2]).queryByText('correct-reference.data.archive.yaml')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectDeactivateEntries[3]).queryByText('lost-reference.data.archive.yaml')).toBeInTheDocument())

  const sectionSelectPaths = screen.getAllByTestId('section-select-path')
  expect(sectionSelectPaths.length).toBe(2)

  await waitFor(() => expect(within(sectionSelectPaths[0]).queryByText('mySubstance1 (./data/mySubstance)')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectPaths[1]).queryByText('mySubstance2 (./data/mySample/mySubstance)')).toBeInTheDocument())

  await act(async () => { userEvent.click(sectionSelectPaths[0]) })
}

const testCreateReferenceDialog = async () => {
  const dialog = screen.getByTestId('create-reference-dialog')
  expect(within(dialog).queryByText('Create new reference')).toBeInTheDocument()

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
    'referenceEditQuantity',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity',
    'AAa34yQsLMIVGUWiQJUxEMK_Fstf',
    'test',
    'password'
  ]
])('test %s', async (name, state, snapshot, entryId, username, password) => {
  await startAPI(state, snapshot, username, password)
  await render(
    <EntryPageContext entryId={entryId}>
      <OverviewView />
    </EntryPageContext>
  )

  await waitForGUI(2000)

  await waitFor(() => expect(screen.getByText('Powder_mixture')).toBeInTheDocument())

  const sectionCards = screen.queryAllByTestId('property-card')
  const powderMixture = sectionCards[0]

  const referenceEditQuantities = within(powderMixture).getAllByTestId('reference-edit-quantity')
  const referenceEditQuantity = referenceEditQuantities[0]
  const inputTextField = within(referenceEditQuantity).getByRole('textbox')
  await waitFor(() => expect(inputTextField.value).toEqual('ref4.archive.yaml#data/mySample/mySubstance'))

  // test section select dialog
  const editReferenceButton = within(referenceEditQuantity).getByTitle('Search for the references').closest('button')
  expect(editReferenceButton).toBeEnabled()

  expect(within(referenceEditQuantity).queryByTitle('Create and assign a new reference')).not.toBeInTheDocument()

  // open edit section select dialog
  await act(async () => { userEvent.click(editReferenceButton) })

  await waitForGUI(1000, true)

  await testSectionSelectDialog()

  // test section select combo
  fireEvent.change(inputTextField, { target: { value: 'ref' } })
  await waitForGUI(2000) // the input changes have been debounced

  await testSectionSelectAutocomplete()

  await waitFor(() => expect(inputTextField.value).toEqual('ref4.archive.yaml#data/mySubstance'))

  // test delete reference
  const deleteReferenceButton = within(referenceEditQuantity).getByTitle('Clear').closest('button')
  expect(deleteReferenceButton).toBeEnabled()
  await act(async () => { await userEvent.click(deleteReferenceButton) })
  await waitFor(() => expect(inputTextField.value).toEqual(''))

  // when the reference is undefined the create reference button should appear
  const createReferenceButton = within(referenceEditQuantity).getByTitle('Create and assign a new reference').closest('button')
  expect(createReferenceButton).toBeEnabled()
  await act(async () => { await userEvent.click(createReferenceButton) })

  await testCreateReferenceDialog()

  await waitForGUI(1000, true)
  await waitFor(() => expect(inputTextField.value).toEqual('newReference.archive.json'))

  closeAPI()
})

test.each([
  [
    'referenceEditQuantity lost entry',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity-lost-entry',
    '6kU3WY1DH5XCpywRJ1lMzNcNL7uX',
    'test',
    'password',
    '',
    'The referenced value does not exist anymore'
  ],
  [
    'referenceEditQuantity lost path',
    'tests.states.entry.references',
    'tests/data/editquantity/referenceEditQuantity-lost-path',
    'GfxzIA3M8X-I719ZuD7wKX_DQJ8S',
    'test',
    'password',
    'ref4.archive.yaml#wrong/path',
    'The provided path does not exist'
  ]
])('test %s', async (name, state, snapshot, entryId, username, password, inputValue, error) => {
  await startAPI(state, snapshot, username, password)
  await render(
    <EntryPageContext entryId={entryId}>
      <OverviewView />
    </EntryPageContext>
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

  within(powderMixture).getByText(error)

  closeAPI()
})
