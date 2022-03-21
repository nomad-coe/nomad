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
  startAPI,
  closeAPI
} from '../conftest'
import { within } from '@testing-library/dom'
import UploadPage from './UploadPage'
import {fireEvent, waitFor} from '@testing-library/react'

const testWritePermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Open the members dialog
  fireEvent.click(screen.getByTestId('edit-members-action'))
  await waitFor(() =>
    expect(screen.queryByText('Main author')).toBeInTheDocument()
  )

  let dialog = screen.getByTestId('edit-members-dialog')
  expect(within(dialog).queryByText('Affiliation')).toBeInTheDocument()
  expect(within(dialog).queryByText('Role')).toBeInTheDocument()

  let rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(3)

  expect(within(rows[0]).queryByText('Markus Scheidgen')).toBeInTheDocument()
  expect(within(rows[0]).queryByText('Main author')).toBeInTheDocument()
  expect(within(rows[0]).getByTestId('member-delete-button')).toBeInTheDocument()
  expect(within(rows[0]).getByTestId('member-delete-button')).toBeDisabled()

  expect(within(rows[1]).queryByText('Sheldon Cooper')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('Co-author')).toBeInTheDocument()
  expect(within(rows[1]).getByTestId('member-delete-button')).toBeInTheDocument()
  expect(within(rows[1]).getByTestId('member-delete-button')).toBeEnabled()

  expect(within(rows[2]).queryByText('Test Tester')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Testeversity')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('Reviewer')).toBeInTheDocument()
  expect(within(rows[2]).getByTestId('member-delete-button')).toBeInTheDocument()
  expect(within(rows[2]).getByTestId('member-delete-button')).toBeEnabled()
}

const testReadPermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Members dialog should be disabled
  expect(screen.getByTestId('edit-members-action')).toBeDisabled()
}

test.each([
  [
    'Published and logged in as main author',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published and logged in as coauthor',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Unpublished and logged in as main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/members-dialog-unpublished',
    'dft_upload',
    'test',
    'password'
  ], [
    'Unpublished and logged in as coauthor',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/members-dialog-unpublished',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Members dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testWritePermissions()
  closeAPI()
})

test.each([
  [
    'Published and logged in as reviewer',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published and logged in as neither reviewer nor writer',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published',
    'dft_upload',
    'admin',
    'password'
  ], [
    'Published and not authenticated',
    'tests.states.uploads.published',
    'tests/data/uploads/members-dialog-published',
    'dft_upload',
    '',
    ''
  ], [
    'Unpublished and logged in as reviewer',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/members-dialog-unpublished',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Members dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testReadPermissions()
  closeAPI()
})
