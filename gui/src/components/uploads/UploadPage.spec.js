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
  waitForGUI,
  startAPI,
  closeAPI
} from '../conftest.spec'
import UploadPage from './UploadPage'
import {fireEvent, waitFor, within, act} from '@testing-library/react'
import userEvent from '@testing-library/user-event'

afterEach(() => closeAPI())

const testShownColumnsAction = async () => {
  await screen.findByTitle('Change the shown columns.')
  expect(screen.queryByText('Parser name')).not.toBeInTheDocument()
  expect(screen.queryByText('Comment')).not.toBeInTheDocument()
  expect(screen.queryByText('References')).not.toBeInTheDocument()
  expect(screen.queryByText('Datasets')).not.toBeInTheDocument()

  const uploadPageRows = screen.queryAllByTestId('datatable-row')
  expect(uploadPageRows.length).toBe(1)
  expect(within(uploadPageRows[0]).queryByText('Mocked')).not.toBeInTheDocument()
  expect(within(uploadPageRows[0]).queryByText('doi')).not.toBeInTheDocument()
  expect(within(uploadPageRows[0]).queryByText('no datasets')).not.toBeInTheDocument()

  act(() => { fireEvent.click(screen.getByTitle('Change the shown columns.')) })
  await waitFor(() => expect(screen.queryByText('Datasets')).toBeInTheDocument())
  await waitFor(() => expect(screen.queryByText('Comment')).toBeInTheDocument())
  await waitFor(() => expect(screen.queryByText('References')).toBeInTheDocument())
  await waitFor(() => expect(screen.queryByText('Parser name')).toBeInTheDocument())

  act(() => { fireEvent.click(screen.getByText('Comment')) })
  act(() => { fireEvent.click(screen.getByText('References')) })
  act(() => { fireEvent.click(screen.getByText('Datasets')) })
  await waitFor(() => expect(within(uploadPageRows[0]).queryByText('Mocked')).toBeInTheDocument())
  await waitFor(() => expect(within(uploadPageRows[0]).queryByText('doi')).toBeInTheDocument())
  await waitFor(() => expect(within(uploadPageRows[0]).queryByText('no datasets')).toBeInTheDocument())
}

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
    'tests/data/uploads/uploadpage-published-author',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published and logged in as coauthor',
    'tests.states.uploads.published',
    'tests/data/uploads/uploadpage-published-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testPublishedWritePermissions()
  await testShownColumnsAction()
})

test.each([
  [
    'Published and logged in as reviewer',
    'tests.states.uploads.published',
    'tests/data/uploads/uploadpage-published-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published and logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published',
    'tests/data/uploads/uploadpage-published-external',
    'dft_upload',
    'admin',
    'password'
  ], [
    'Published and not authenticated',
    'tests.states.uploads.published',
    'tests/data/uploads/uploadpage-published-noauth',
    'dft_upload',
    '',
    ''
  ], [
    'Unpublished and logged in as reviewer',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/uploadpage-unpublished-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published with embargo and logged in as reviewer',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/uploadpage-published-with-embargo-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testReadOnlyPermissions()
  await testShownColumnsAction()
})

test.each([
  [
    'Unpublished and logged in as main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/uploadpage-unpublished-author',
    'dft_upload',
    'test',
    'password'
  ], [
    'Unpublished and logged in as coauthor',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/uploadpage-unpublished-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testUnpublishedWritePermissions()
  await testShownColumnsAction()
})

test('Render upload page: multiple entries', async () => {
  await startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/multiple_entries', 'test', 'password')
  render(<UploadPage uploadId="dft_upload_1"/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table header is rendered correctly
  expect(screen.queryByText('6 entries')).toBeInTheDocument()
  expect(screen.queryByTestId('table-pagination')).toBeInTheDocument()
  expect(screen.queryByTestId('datatable-body')).toBeInTheDocument()

  const datatableBody = screen.getByTestId('datatable-body')

  // Test if the pagination works correctly
  const rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)
  expect(within(datatableBody).queryByText('vasp_6.xml')).not.toBeInTheDocument()

  // Test if the name of the entries are rendered in the right order
  expect(within(rows[0]).queryByText('vasp_1.xml')).toBeInTheDocument()
  expect(within(rows[1]).queryByText('vasp_2.xml')).toBeInTheDocument()
  expect(within(rows[2]).queryByText('vasp_3.xml')).toBeInTheDocument()
  expect(within(rows[3]).queryByText('vasp_4.xml')).toBeInTheDocument()
  expect(within(rows[4]).queryByText('vasp_5.xml')).toBeInTheDocument()
})

test('Delete selected entries from table', async () => {
  await startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/delete_entries_from_table', 'test', 'password')
  render(<UploadPage uploadId="dft_upload_1"/>)

  // Wait for the page to load, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table header is rendered correctly
  expect(screen.queryByText('6 entries')).toBeInTheDocument()
  expect(screen.queryByTestId('table-pagination')).toBeInTheDocument()
  let rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(5)

  // TODO: temporary workaround, the page updates several times for some reasons, and if we're
  // clicking too quickly (before the initial cascade of rendering is completely done), this can
  // result in problems.
  await waitForGUI()

  // Go to the next page in the entry table
  const nextPage = screen.getByButtonText(/next page/i)
  await userEvent.click(nextPage)
  await waitFor(() => {
    const table = screen.getByTestId('datatable-body')
    expect(within(table).getByText('vasp_6.xml')).toBeVisible()
    expect(within(table).queryAllByText('vasp_1.xml').length).toBe(0)
  })

  expect(screen.queryAllByButtonText('Delete selected entries').length).toBe(0)

  // Select the last entry
  rows = screen.queryAllByTestId('datatable-row')
  await userEvent.click(within(rows[0]).getByRole('checkbox'))

  // Wait for delete entries button to appear
  let deleteButton = await screen.findByButtonText('Delete selected entries')

  // Delete the entry
  await userEvent.click(deleteButton)
  let deleteConfirmButton = await screen.findByButtonText('Delete mainfiles')
  await userEvent.click(deleteConfirmButton)
  // Should delete and go back to the first page
  await waitFor(() => {
    expect(screen.queryByText('5 entries')).toBeInTheDocument()
    expect(screen.queryAllByTestId('datatable-row').length).toBe(5)
    const table = screen.getByTestId('datatable-body')
    expect(within(table).getByText('vasp_1.xml')).toBeVisible()
    expect(within(table).queryAllByText('vasp_6.xml').length).toBe(0)
  })

  // Select two entries (#3 and #5)
  rows = screen.queryAllByTestId('datatable-row')
  await userEvent.click(within(rows[2]).getByRole('checkbox'))
  await userEvent.click(within(rows[4]).getByRole('checkbox'))

  // Wait for delete entries button to appear
  deleteButton = await screen.findByButtonText('Delete selected entries')

  // Delete the entries
  await userEvent.click(deleteButton)
  deleteConfirmButton = await screen.findByButtonText('Delete mainfiles')
  await userEvent.click(deleteConfirmButton)
  await waitFor(() => {
    expect(screen.queryByText('3 entries')).toBeInTheDocument()
    expect(screen.queryAllByTestId('datatable-row').length).toBe(3)
    const table = screen.getByTestId('datatable-body')
    expect(within(table).getByText('vasp_1.xml')).toBeVisible()
    expect(within(table).getByText('vasp_2.xml')).toBeVisible()
    expect(within(table).getByText('vasp_4.xml')).toBeVisible()
    expect(within(table).queryAllByText('vasp_3.xml').length).toBe(0)
    expect(within(table).queryAllByText('vasp_5.xml').length).toBe(0)
  })
})

test.each([
  [
    'Published with embargo and logged in as main author',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/uploadpage-published-with-embargo-author',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published with embargo and logged in as coauthor',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/uploadpage-published-with-embargo-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testEmbargoedPublishesWritePermissions()
  await testShownColumnsAction()
})

test.each([
  [
    'unpublished, not authenticated',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/uploadpage-unpublished-noauth',
    'dft_upload',
    '',
    '',
    'You need to log in to access the specified upload.'
  ],
  [
    'unpublished, logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/uploadpage-unpublished-external',
    'dft_upload',
    'admin',
    'password',
    'You do not have access to the specified upload.'
  ],
  [
    'published with embargo, not authenticated',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/uploadpage-published-with-embargo-noauth',
    'dft_upload',
    '',
    '',
    'You do not have access to the specified upload - published with embargo.'
  ], [
    'published with embargo, logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/uploadpage-published-with-embargo-external',
    'dft_upload',
    'admin',
    'password',
    'You do not have access to the specified upload - published with embargo.'
  ], [
    'unknown upload_id',
    'tests.states.uploads.published',
    'tests/data/uploads/uploadpage-nonexistent',
    'a_not_exist_upload_ID',
    '',
    '',
    'The specified upload_id was not found.'
  ]
])('Render upload page: error message due to %s', async (name, state, snapshot, uploadId, username, password, msg) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await screen.findByText(msg)
})
