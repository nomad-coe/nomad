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

import { waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import React from 'react'
import { closeAPI, render, screen, startAPI, within } from '../conftest.spec'
import LoginLogout from '../LoginLogout'
import UploadPage from './UploadPage'

afterEach(() => {
  closeAPI()
})

// MUI's Select component uses divs etc. and mounts the options temporarily near the root
const selectRole = async (row, initialRole, targetRole, user) => {
  const roleInput = within(row).getByDisplayValue(initialRole)

  const roleButton = within(row).getByRoleAndText('button', initialRole)
  await user.click(roleButton)

  const roleOptionList = await screen.getByRoleAndText('listbox', initialRole)
  const roleTargetOption = within(roleOptionList).getByText(targetRole)
  await user.click(roleTargetOption)

  expect(roleInput.value).toBe(targetRole)
}

const searchAndAddMember = async (autocompleteInput, searchText, memberText, user) => {
  await user.click(autocompleteInput)
  await screen.findByText('No options')
  await user.type(autocompleteInput, searchText)
  expect(autocompleteInput.value).toEqual(searchText)
  await screen.findByText('Loading…')
  if (memberText) {
    await screen.findByRoleAndText('option', memberText)
    await user.keyboard('[ArrowDown]')
    await user.keyboard('[Enter]')
  } else {
    await screen.findByText('No options')
    await user.tab()
  }
  expect(autocompleteInput.value).toEqual('')
}

const checkRow = (row, expectedValues) => {
  const { texts, isDeletable = true } = expectedValues
  texts.forEach((text) => expect(within(row).getByText(text)).toBeInTheDocument())
  const deleteButton = within(row).getByTestId('member-delete-button')
  isDeletable ? expect(deleteButton).toBeEnabled() : expect(deleteButton).toBeDisabled()
}

const testInitialDialog = (dialog) => {
  const rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)
  checkRow(rows[0], { texts: ['Markus Scheidgen', 'Main author'], isDeletable: false })
  checkRow(rows[1], { texts: ['Sheldon Cooper', 'Testeversity', 'Co-author']})
  checkRow(rows[2], { texts: ['Test Tester', 'Testeversity', 'Reviewer']})
  checkRow(rows[3], { texts: ['Group Cooper', 'Sheldon Cooper', 'Co-author']})
  checkRow(rows[4], { texts: ['Group Tester', 'Test Tester', 'Reviewer']})
}

const testNewDialog = (dialog) => {
  const rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)
  checkRow(rows[0], { texts: ['Markus Scheidgen', 'Main author'], isDeletable: false })
  checkRow(rows[1], { texts: ['Admin Administrator', 'Co-author']})
  checkRow(rows[2], { texts: ['Test Tester', 'Testeversity', 'Reviewer']})
  checkRow(rows[3], { texts: ['Group Admin', 'Admin Administrator', 'Co-author']})
  checkRow(rows[4], { texts: ['Group Tester', 'Test Tester', 'Reviewer']})
}

const testAddRemoveMembers = async (dialog, user) => {
  const searchMembers = within(dialog).getByRole('combobox')
  const autocompleteInput = within(searchMembers).getByRole('textbox')

  // search user that is already in the list, no options should be shown
  await searchAndAddMember(autocompleteInput, 'teste', null, user)
  expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(5)

  // search user that is not in the list and add it
  await searchAndAddMember(autocompleteInput, 'admin', 'Admin Administrator', user)
  await waitFor(async () => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(6))

  // search group that is not in the list and add it
  const searchTypeSelect = within(dialog).getByRole('button', {'name': 'User'})
  await user.click(searchTypeSelect)
  await user.click(screen.getByText('Group'))
  expect(searchTypeSelect).toHaveTextContent('Group')

  await searchAndAddMember(autocompleteInput, 'admin', 'Admin Administrator', user)
  await waitFor(async () => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(7))

  let rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(7)
  checkRow(rows[0], { texts: ['Markus Scheidgen', 'Main author'], isDeletable: false })
  checkRow(rows[1], { texts: ['Sheldon Cooper', 'Testeversity', 'Co-author']})
  checkRow(rows[2], { texts: ['Test Tester', 'Testeversity', 'Reviewer']})
  checkRow(rows[3], { texts: ['Group Cooper', 'Sheldon Cooper', 'Co-author']})
  checkRow(rows[4], { texts: ['Group Tester', 'Test Tester', 'Reviewer']})
  checkRow(rows[5], { texts: ['Admin Administrator', 'Reviewer']})
  checkRow(rows[6], { texts: ['Group Admin', 'Admin Administrator', 'Reviewer']})

  await user.click(within(rows[1]).getByTestId('member-delete-button'))
  await waitFor(async () => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(6))

  rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(6)
  checkRow(rows[0], { texts: ['Markus Scheidgen', 'Main author'], isDeletable: false })
  checkRow(rows[1], { texts: ['Test Tester', 'Testeversity', 'Reviewer']})
  checkRow(rows[2], { texts: ['Group Cooper', 'Sheldon Cooper', 'Co-author']})
  checkRow(rows[3], { texts: ['Group Tester', 'Test Tester', 'Reviewer']})
  checkRow(rows[4], { texts: ['Admin Administrator', 'Reviewer']})
  checkRow(rows[5], { texts: ['Group Admin', 'Admin Administrator', 'Reviewer']})

  await user.click(within(rows[2]).getByTestId('member-delete-button'))
  await waitFor(async () => expect(within(dialog).queryAllByTestId('datatable-row').length).toBe(5))

  rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)
  checkRow(rows[0], { texts: ['Markus Scheidgen', 'Main author'], isDeletable: false })
  checkRow(rows[1], { texts: ['Test Tester', 'Testeversity', 'Reviewer']})
  checkRow(rows[2], { texts: ['Group Tester', 'Test Tester', 'Reviewer']})
  checkRow(rows[3], { texts: ['Admin Administrator', 'Reviewer']})
  checkRow(rows[4], { texts: ['Group Admin', 'Admin Administrator', 'Reviewer']})

  await selectRole(rows[3], 'Reviewer', 'Co-author', user)
  await selectRole(rows[4], 'Reviewer', 'Co-author', user)

  rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)
  checkRow(rows[0], { texts: ['Markus Scheidgen', 'Main author'], isDeletable: false })
  checkRow(rows[1], { texts: ['Test Tester', 'Testeversity', 'Reviewer']})
  checkRow(rows[2], { texts: ['Group Tester', 'Test Tester', 'Reviewer']})
  checkRow(rows[3], { texts: ['Admin Administrator', 'Co-author']})
  checkRow(rows[4], { texts: ['Group Admin', 'Admin Administrator', 'Co-author']})
}

const submitChanges = async (dialog, user) => {
  const submitButton = within(dialog).queryByText('Submit')
  expect(submitButton).toBeEnabled()
  await user.click(submitButton)
  await waitFor(() => expect(dialog).not.toBeInTheDocument())
  const processing = await screen.findByText('Upload is processing ...')
  await waitFor(() => expect(processing).not.toBeInTheDocument())
}

const openMembersDialog = async (user) => {
  await screen.findByTestId('logout-button')
  await screen.findByText('unnamed upload')
  const editMembersButton = screen.getByTestId('edit-members-action')
  await waitFor(() => expect(editMembersButton).toBeEnabled())
  await user.click(editMembersButton)
  const dialog = await screen.findByTestId('edit-members-dialog')
  expect(within(dialog).getByText('Name')).toBeInTheDocument()
  expect(within(dialog).getByText('Affiliation')).toBeInTheDocument()
  expect(within(dialog).getByText('Role')).toBeInTheDocument()
  await within(dialog).findByText('Main author')
  return dialog
}

const testReadOnlyPermissions = async (isLoggedIn) => {
  await screen.findByTestId(isLoggedIn ? 'logout-button' : 'login-register-button')
  await screen.findByText('unnamed upload')
  expect(screen.getByTestId('edit-members-action')).toBeDisabled()
}

test.each([
  [
    'Published and logged in as main author',
    'tests.states.uploads.published_twin_access',
    'tests/data/uploads/members-dialog-published-author',
    'dft_upload',
    'test',
    'password'
  ],
  [
    'Published and logged in as coauthor',
    'tests.states.uploads.published_twin_access',
    'tests/data/uploads/members-dialog-published-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ],
  [
    'Unpublished and logged in as main author',
    'tests.states.uploads.unpublished_twin_access',
    'tests/data/uploads/members-dialog-unpublished-author',
    'dft_upload',
    'test',
    'password'
  ]
])('Members dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  const user = userEvent.setup()
  render(<><LoginLogout/><UploadPage uploadId={uploadId}/></>)

  const dialog = await openMembersDialog(user)
  testInitialDialog(dialog)
  await testAddRemoveMembers(dialog, user)
  await submitChanges(dialog, user)

  if (username === 'test') {
    const newDialog = await openMembersDialog(user)
    await testNewDialog(newDialog)
  } else if (username === 'scooper') {
    await testReadOnlyPermissions(true)
  }
})

test.each([
  [
    'Published and logged in as reviewer',
    'tests.states.uploads.published_twin_access',
    'tests/data/uploads/members-dialog-published-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published and logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published_twin_access',
    'tests/data/uploads/members-dialog-published-external',
    'dft_upload',
    'admin',
    'password'
  ], [
    'Published and not authenticated',
    'tests.states.uploads.published_twin_access',
    'tests/data/uploads/members-dialog-published-noauth',
    'dft_upload',
    '',
    ''
  ], [
    'Unpublished and logged in as reviewer',
    'tests.states.uploads.unpublished_twin_access',
    'tests/data/uploads/members-dialog-unpublished-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Members dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<><LoginLogout/><UploadPage uploadId={uploadId}/></>)
  const isLoggedIn = username !== ''
  await testReadOnlyPermissions(isLoggedIn)
})
