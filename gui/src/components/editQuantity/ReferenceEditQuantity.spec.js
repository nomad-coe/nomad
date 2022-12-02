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
import {render, screen, startAPI, closeAPI} from '../conftest.spec'
import OverviewView from '../entry/OverviewView'
import EntryPageContext from '../entry/EntryPageContext'
import { act } from 'react-dom/test-utils'
import userEvent from '@testing-library/user-event'
import {fireEvent} from '@testing-library/react'

const testSectionSelectDialog = async () => {
  const dialog = screen.getByTestId('section-select-dialog')
  expect(within(dialog).queryByText('Filters')).toBeInTheDocument()

  await waitFor(() => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(4))
  const rows = within(dialog).getAllByTestId('datatable-row')

  await waitFor(() => expect(within(rows[0]).queryByText('Copper (II) Selenide')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[0]).queryByText('Chemical')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[1]).queryByText('./data')).toBeInTheDocument())

  const sections = within(rows[1]).queryAllByRole('menuitem')
  expect(sections.length).toBe(1)
  await waitFor(() => expect(within(sections[0]).queryByText('./data')).toBeInTheDocument())
  await waitFor(() => expect(within(sections[0]).getByTestId('check-icon')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[2]).queryByText('Tin (II) Selenide')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[2]).queryByText('Chemical')).toBeInTheDocument())

  await waitFor(() => expect(within(rows[3]).queryByText('Zinc Selenide')).toBeInTheDocument())
  await waitFor(() => expect(within(rows[3]).queryByText('Chemical')).toBeInTheDocument())

  // close the dialog
  await userEvent.click(within(dialog).getByRole('button', {name: /cancel/i}))
  await waitFor(() => expect(screen.queryByTestId('section-select-dialog')).not.toBeInTheDocument())
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
    'test referenceEditQuantity',
    'tests.states.entry.eln',
    'tests/data/editquantity/referenceEditQuantity',
    'bC7byHvWJp62Sn9uiuJUB38MT5j-',
    'test',
    'password'
  ]
])('eln %s', async (name, state, snapshot, entryId, username, password) => {
  await startAPI(state, snapshot, username, password)
  await act(async () => render(
    <EntryPageContext entryId={entryId}>
      <OverviewView />
    </EntryPageContext>
  ))

  await waitFor(() => expect(screen.getByText('HotplateAnnealing')).toBeInTheDocument())

  const sectionCards = screen.queryAllByTestId('property-card')
  const cardSample = sectionCards[0]
  const cardPvdEvaporation = sectionCards[1]

  const referenceEditQuantities = within(cardSample).getAllByTestId('reference-edit-quantity')
  const referenceEditQuantity = referenceEditQuantities[0]
  const inputTextField = within(referenceEditQuantity).getByRole('textbox')
  await waitFor(() => expect(inputTextField.value).toEqual('Copper_II_Selenide.archive.json'))

  const editReferenceButton = within(referenceEditQuantity).getByTitle('Search for the references').closest('button')
  expect(editReferenceButton).toBeEnabled()

  expect(within(referenceEditQuantity).queryByTitle('Create and assign a new reference')).not.toBeInTheDocument()

  // open edit section select dialog
  await userEvent.click(editReferenceButton)

  await testSectionSelectDialog()

  // test delete reference
  const pvdReferenceEditQuantities = within(cardPvdEvaporation).getAllByTestId('reference-edit-quantity')
  const pvdReferenceEditQuantity = pvdReferenceEditQuantities[0]
  const pvdInputTextField = within(pvdReferenceEditQuantity).getByRole('textbox')
  await waitFor(() => expect(pvdInputTextField.value).toEqual('PVD-P.archive.json'))

  const deleteReferenceButton = within(pvdReferenceEditQuantity).getByTitle('Delete the reference').closest('button')
  expect(deleteReferenceButton).toBeEnabled()
  await userEvent.click(deleteReferenceButton)
  await waitFor(() => expect(pvdInputTextField.value).toEqual(''))

  // when the reference is undefined the create reference button should appear
  const createReferenceButton = within(pvdReferenceEditQuantity).getByTitle('Create and assign a new reference').closest('button')
  expect(createReferenceButton).toBeEnabled()
  await userEvent.click(createReferenceButton)

  await testCreateReferenceDialog()

  await waitFor(() => expect(pvdInputTextField.value).toEqual('newReference.archive.json'))

  closeAPI()
})
