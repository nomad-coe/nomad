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
} from '../conftest.spec'
import { within } from '@testing-library/dom'
import UploadPage from './UploadPage'

const testPublishedWritePermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if only the first two steps are shown
  expect(screen.queryByText('Prepare and upload your files')).toBeInTheDocument()
  expect(screen.queryByText('Processing completed, 1/1 entries processed')).toBeInTheDocument()
  expect(screen.queryByText('You can either select and edit individual entries from the list above, or edit all entries at once.')).toBeInTheDocument()
  expect(screen.queryByText('This upload has already been published.')).toBeInTheDocument()

  expect(screen.getByTestId('edit-members-action')).toBeEnabled()
  expect(screen.getByTestId('upload-download-action')).toBeEnabled()
  expect(screen.getByTestId('upload-reprocess-action')).toBeDisabled()
  expect(screen.getByTestId('source-api-action')).toBeEnabled()
  expect(screen.getByTestId('upload-delete-action')).toBeDisabled()
}

const testUnpublishedWritePermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table title is rendered correctly
  expect(screen.queryByText('1 entry')).toBeInTheDocument()

  // Test if only the first two steps are shown
  expect(screen.queryByText('Prepare and upload your files')).toBeInTheDocument()
  expect(screen.queryByText('Processing completed, 1/1 entries processed')).toBeInTheDocument()
  expect(screen.queryByText('You can either select and edit individual entries from the list above, or edit all entries at once.')).toBeInTheDocument()

  expect(screen.getByTestId('edit-members-action')).toBeEnabled()
  expect(screen.getByTestId('upload-download-action')).toBeEnabled()
  expect(screen.getByTestId('upload-reprocess-action')).toBeEnabled()
  expect(screen.getByTestId('source-api-action')).toBeEnabled()
  expect(screen.getByTestId('upload-delete-action')).toBeEnabled()
  expect(screen.getByTestId('edit-metadata-button')).toBeEnabled()
  expect(screen.getByTestId('publish-upload-button')).toBeEnabled()
}

const testEmbargoedPublishesWritePermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if only the first two steps are shown
  expect(screen.queryByText('Prepare and upload your files')).toBeInTheDocument()
  expect(screen.queryByText('Processing completed, 1/1 entries processed')).toBeInTheDocument()
  expect(screen.queryByText('You can either select and edit individual entries from the list above, or edit all entries at once.')).toBeInTheDocument()

  expect(screen.getByTestId('edit-members-action')).toBeEnabled()
  expect(screen.getByTestId('upload-download-action')).toBeEnabled()
  expect(screen.getByTestId('upload-reprocess-action')).toBeDisabled()
  expect(screen.getByTestId('source-api-action')).toBeEnabled()
  expect(screen.getByTestId('upload-delete-action')).toBeDisabled()
  expect(screen.getByTestId('edit-metadata-button')).toBeEnabled()
  expect(screen.getByText('Lift Embargo')).toBeEnabled()
}

const testReadOnlyPermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table title is rendered correctly
  expect(screen.queryByText('1 entry')).toBeInTheDocument()

  // Test if only the first two steps are shown
  expect(screen.queryByText('Prepare and upload your files')).toBeInTheDocument()
  expect(screen.queryByText('Processing completed, 1/1 entries processed')).toBeInTheDocument()
  expect(screen.queryByText('You can either select and edit individual entries from the list above, or edit all entries at once.')).not.toBeInTheDocument()

  expect(screen.queryByTestId('edit-members-action')).toBeDisabled()
  expect(screen.queryByTestId('upload-download-action')).toBeEnabled()
  expect(screen.queryByTestId('upload-reprocess-action')).toBeDisabled()
  expect(screen.queryByTestId('source-api-action')).toBeEnabled()
  expect(screen.queryByTestId('upload-delete-action')).toBeDisabled()
  expect(screen.queryByTestId('edit-metadata-button')).not.toBeInTheDocument()
  expect(screen.queryByTestId('publish-upload-button')).not.toBeInTheDocument()
}

test.each([
  [
    'Published and logged in as main author',
    'tests.states.uploads.published',
    'tests/data/uploads/published',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published and logged in as coauthor',
    'tests.states.uploads.published',
    'tests/data/uploads/published',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testPublishedWritePermissions()
  closeAPI()
})

test.each([
  [
    'Published and logged in as reviewer',
    'tests.states.uploads.published',
    'tests/data/uploads/published',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published and logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published',
    'tests/data/uploads/published',
    'dft_upload',
    'admin',
    'password'
  ], [
    'Published and not authenticated',
    'tests.states.uploads.published',
    'tests/data/uploads/published',
    'dft_upload',
    '',
    ''
  ], [
    'Unpublished and logged in as reviewer',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/unpublished',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published with embargo and logged in as reviewer',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/published_with_embargo',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testReadOnlyPermissions()
  closeAPI()
})

test.each([
  [
    'Unpublished and logged in as main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/unpublished',
    'dft_upload',
    'test',
    'password'
  ], [
    'Unpublished and logged in as coauthor',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/unpublished',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testUnpublishedWritePermissions()
  closeAPI()
})

test('Render upload page: multiple entries', async () => {
  startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/multiple_entries')
  render(<UploadPage uploadId="dft_upload_1"/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table header is rendered correctly
  expect(screen.queryByText('6 entries')).toBeInTheDocument()
  expect(screen.queryByTestId('table-pagination')).toBeInTheDocument()
  expect(screen.queryByTestId('datatable-body')).toBeInTheDocument()

  let datatableBody = screen.getByTestId('datatable-body')

  // Test if the pagination works correctly
  let rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)
  expect(within(datatableBody).queryByText('vasp_6.xml')).not.toBeInTheDocument()

  // Test if the name of the entries are rendered in the right order
  expect(within(rows[0]).queryByText('vasp_1.xml')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('vasp_2.xml')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('vasp_3.xml')).toBeInTheDocument()
  expect(within(rows[3]).queryByText('vasp_4.xml')).toBeInTheDocument()
  expect(within(rows[4]).queryByText('vasp_5.xml')).toBeInTheDocument()

  closeAPI()
})

test.each([
  [
    'Published with embargo and logged in as main author',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/published_with_embargo',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published with embargo and logged in as coauthor',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/published_with_embargo',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testEmbargoedPublishesWritePermissions()
  closeAPI()
})

test.each([
  [
    'unpublished, not authenticated',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/unpublished-not-writer',
    'dft_upload',
    '',
    '',
    'You do not have access to the specified upload - not published yet.'
  ],
  [
    'unpublished, logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/unpublished-not-writer',
    'dft_upload',
    'admin',
    'password',
    'You do not have access to the specified upload - not published yet.'
  ],
  [
    'published with embargo, not authenticated',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/published_with_embargo-not-writer',
    'dft_upload',
    '',
    '',
    'You do not have access to the specified upload - published with embargo.'
  ], [
    'published with embargo, logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/published_with_embargo-not-writer',
    'dft_upload',
    'admin',
    'password',
    'You do not have access to the specified upload - published with embargo.'
  ], [
    'unknown upload_id',
    'tests.states.uploads.published',
    'tests/data/uploads/not_exists',
    'a_not_exist_upload_ID',
    '',
    '',
    'The specified upload_id was not found.'
  ]
])('Render upload page: error message due to %s', async (name, state, snapshot, uploadId, username, password, msg) => {
  startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await screen.findByText(msg)
  closeAPI()
})
