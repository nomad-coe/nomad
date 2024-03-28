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
import {
  render,
  screen,
  within,
  startAPI,
  closeAPI, waitForGUI
} from '../conftest.spec'
import UploadPage from './UploadPage'
import {act, fireEvent, waitFor} from '@testing-library/react'

const testInitialDialog = async (dialog) => {
  const rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(3)

  expect(within(rows[0]).queryByText('Markus Scheidgen')).toBeInTheDocument()
  expect(within(rows[0]).queryByText('Main author')).toBeInTheDocument()
  expect(within(rows[0]).getByTestId('member-delete-button')).toBeDisabled()

  expect(within(rows[1]).queryByText('Sheldon Cooper')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Co-author')).toBeInTheDocument()
  expect(within(rows[1]).getByTestId('member-delete-button')).toBeEnabled()

  expect(within(rows[2]).queryByText('Test Tester')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Reviewer')).toBeInTheDocument()
  expect(within(rows[2]).getByTestId('member-delete-button')).toBeEnabled()
}

const testNewDialog = async (dialog) => {
  const rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(3)

  expect(within(rows[0]).queryByText('Markus Scheidgen')).toBeInTheDocument()
  expect(within(rows[0]).queryByText('Main author')).toBeInTheDocument()
  expect(within(rows[0]).getByTestId('member-delete-button')).toBeDisabled()

  expect(within(rows[1]).queryByText('Admin Administrator')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Co-author')).toBeInTheDocument()
  expect(within(rows[1]).getByTestId('member-delete-button')).toBeEnabled()

  expect(within(rows[2]).queryByText('Test Tester')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Reviewer')).toBeInTheDocument()
  expect(within(rows[2]).getByTestId('member-delete-button')).toBeEnabled()
}

const testAddRemoveMembers = async (dialog) => {
  const searchMembers = within(dialog).queryAllByRole('combobox')[0]
  const autocompleteInput = within(searchMembers).getByRole('textbox')
  const addMemberButton = within(dialog).getByButtonText('Add')

  await waitFor(() => expect(within(dialog).queryByText('The selected user is already in the members list')).not.toBeVisible())

  searchMembers.focus()
  // assign an incomplete value to the input field
  fireEvent.change(autocompleteInput, { target: { value: 'teste' } })
  await waitForGUI(700, true)
  await waitFor(() => expect(autocompleteInput.value).toEqual('teste'))
  await waitForGUI(3500, true)
  fireEvent.keyDown(searchMembers, { key: 'ArrowDown' })
  fireEvent.keyDown(searchMembers, { key: 'Enter' })
  await waitForGUI()
  // test if it is completed
  await waitFor(() => expect(autocompleteInput.value).toEqual('Test Tester (Testeversity)'))
  await waitFor(() => expect(addMemberButton).toBeDisabled())
  await waitFor(() => expect(within(dialog).queryByText('The selected user is already in the members list')).toBeVisible())

  searchMembers.focus()
  // assign an incomplete value to the input field
  fireEvent.change(autocompleteInput, { target: { value: '' } })
  fireEvent.change(autocompleteInput, { target: { value: 'admin' } })
  await waitForGUI(700, true)
  await waitFor(() => expect(autocompleteInput.value).toEqual('admin'))
  await waitForGUI(3500, true)
  fireEvent.keyDown(searchMembers, { key: 'ArrowDown' })
  fireEvent.keyDown(searchMembers, { key: 'Enter' })
  await waitForGUI()
  await waitFor(() => expect(autocompleteInput.value).toEqual('Admin Administrator'))
  await waitFor(() => expect(addMemberButton).toBeEnabled())
  await waitFor(() => expect(within(dialog).queryByText('The selected user is already in the members list')).not.toBeVisible())

  fireEvent.click(addMemberButton)
  await waitFor(() => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(4))

  const rows = within(dialog).queryAllByTestId('datatable-row')
  expect(within(rows[0]).queryByText('Markus Scheidgen')).toBeInTheDocument()
  expect(within(rows[0]).queryByText('Main author')).toBeInTheDocument()
  expect(within(rows[0]).getByTestId('member-delete-button')).toBeDisabled()

  expect(within(rows[1]).queryByText('Sheldon Cooper')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Co-author')).toBeInTheDocument()
  expect(within(rows[1]).getByTestId('member-delete-button')).toBeEnabled()

  expect(within(rows[2]).queryByText('Test Tester')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Reviewer')).toBeInTheDocument()
  expect(within(rows[2]).getByTestId('member-delete-button')).toBeEnabled()

  expect(within(rows[3]).queryByText('Admin Administrator')).toBeInTheDocument()
  expect(within(rows[3]).queryByText('Co-author')).toBeInTheDocument()
  expect(within(rows[3]).getByTestId('member-delete-button')).toBeEnabled()

  fireEvent.click(within(rows[1]).getByTestId('member-delete-button'))
  await waitFor(() => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(3))
}

const submitChanges = async (dialog) => {
  const submitButton = within(dialog).queryByText('Submit')
  expect(submitButton).toBeEnabled()
  fireEvent.click(submitButton)
  await waitFor(() => expect(screen.queryByTestId('edit-members-dialog')).not.toBeInTheDocument())
}

const openMembersDialog = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Open the members dialog
  await act(async () => { fireEvent.click(screen.getByTestId('edit-members-action')) })
  await waitFor(() => expect(screen.queryByText('Main author')).toBeInTheDocument())
  const dialog = screen.getByTestId('edit-members-dialog')
  expect(within(dialog).queryByText('Affiliation')).toBeInTheDocument()
  expect(within(dialog).queryByText('Role')).toBeInTheDocument()
  return dialog
}

const testReadOnlyPermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Members dialog should be disabled
  expect(screen.getByTestId('edit-members-action')).toBeDisabled()
}

test.each([
  [
    'Published and logged in as main author',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published-author',
    'dft_upload',
    'test',
    'password'
  ],
  [
    'Published and logged in as coauthor',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ],
  [
    'Unpublished and logged in as main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/members-dialog-unpublished-author',
    'dft_upload',
    'test',
    'password'
  ]
])('Members dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)

  const dialog = await openMembersDialog()
  await testInitialDialog(dialog)
  await testAddRemoveMembers(dialog)

  await submitChanges(dialog)
  await waitForGUI(2000, true)

  if (username === 'test') {
    const newDialog = await openMembersDialog()
    waitForGUI()
    await testNewDialog(newDialog)
  } else if (username === 'scooper') {
    await testReadOnlyPermissions()
  }

  closeAPI()
})

test.each([
  [
    'Published and logged in as reviewer',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published and logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published-external',
    'dft_upload',
    'admin',
    'password'
  ], [
    'Published and not authenticated',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published-noauth',
    'dft_upload',
    '',
    ''
  ], [
    'Unpublished and logged in as reviewer',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/members-dialog-unpublished-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Members dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testReadOnlyPermissions()
  closeAPI()
})

test('Toggle visible for all checkbox', async () => {
  await startAPI(
    'tests.states.uploads.unpublished',
    'tests/data/uploads/members-dialog-toggle-visible',
    'test',
    'password'
  )
  waitForGUI(0, true)
  render(<UploadPage uploadId='dft_upload'/>)

  let dialog = await openMembersDialog()
  let checkbox = await within(dialog).findByRoleAndText('checkbox', 'Visible for all')
  expect(checkbox.checked).toEqual(false)
  fireEvent.click(checkbox)
  expect(checkbox.checked).toEqual(true)
  await submitChanges(dialog)
  await waitForGUI(2000, true)

  dialog = await openMembersDialog()
  checkbox = await within(dialog).findByRoleAndText('checkbox', 'Visible for all')
  expect(checkbox.checked).toEqual(true)
  fireEvent.click(checkbox)
  expect(checkbox.checked).toEqual(false)
  await submitChanges(dialog)
  await waitForGUI(2000, true)

  dialog = await openMembersDialog()
  checkbox = await within(dialog).findByRoleAndText('checkbox', 'Visible for all')
  expect(checkbox.checked).toEqual(false)

  closeAPI()
})
