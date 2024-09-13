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

import { act, fireEvent, waitFor, within } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import React from 'react'
import { closeAPI, render, screen, startAPI, waitForGUI } from '../conftest.spec'
import LoginLogout from '../LoginLogout'
import UploadPage from './UploadPage'

afterEach(() => closeAPI())

const queryVisibleForAllCheckbox = () => {
  return screen.queryByRole('checkbox', { name: /Enabling this will allow all users/ })
}

const testShownColumnsAction = async () => {
  await screen.findByTitle('Change the shown columns')
  const processingTable = screen.getByTestId('processing-table')
  expect(within(processingTable).queryByText('Parser name')).not.toBeInTheDocument()
  expect(within(processingTable).queryByText('Comment')).not.toBeInTheDocument()
  expect(within(processingTable).queryByText('References')).not.toBeInTheDocument()
  expect(within(processingTable).queryByText('Datasets')).not.toBeInTheDocument()

  const uploadPageRows = screen.queryAllByTestId('datatable-row')
  expect(uploadPageRows.length).toBe(1)
  expect(within(uploadPageRows[0]).queryByText('Mocked')).not.toBeInTheDocument()
  expect(within(uploadPageRows[0]).queryByText('doi')).not.toBeInTheDocument()
  expect(within(uploadPageRows[0]).queryByText('no datasets')).not.toBeInTheDocument()

  act(() => { fireEvent.click(screen.getByTitle('Change the shown columns')) })
  const selectMenu = screen.getByTestId('column-select-menu')
  await waitFor(() => expect(within(selectMenu).queryByText('Datasets')).toBeInTheDocument())
  await waitFor(() => expect(within(selectMenu).queryByText('Comment')).toBeInTheDocument())
  await waitFor(() => expect(within(selectMenu).queryByText('References')).toBeInTheDocument())
  await waitFor(() => expect(within(selectMenu).queryByText('Parser name')).toBeInTheDocument())

  act(() => { fireEvent.click(within(selectMenu).getByText('Comment')) })
  act(() => { fireEvent.click(within(selectMenu).getByText('References')) })
  act(() => { fireEvent.click(within(selectMenu).getByText('Datasets')) })
  await waitFor(() => expect(within(uploadPageRows[0]).queryByText('Mocked')).toBeInTheDocument())
  await waitFor(() => expect(within(uploadPageRows[0]).queryByText('doi')).toBeInTheDocument())
  await waitFor(() => expect(within(uploadPageRows[0]).queryByText('no datasets')).toBeInTheDocument())
}

const testPublishedWritePermissions = async () => {
  await screen.findByTestId('logout-button')
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if only the first two steps are shown
  expect(screen.queryByText('Prepare and upload your files')).toBeInTheDocument()
  expect(screen.queryByText('Processing completed, 1/1 entries processed')).toBeInTheDocument()
  expect(screen.queryByText('You can either select and edit individual entries from the list above, or edit all entries at once.')).toBeInTheDocument()
  expect(screen.queryByText('This upload has already been published.')).toBeInTheDocument()

  expect(screen.getByTestId('edit-members-action')).toBeEnabled()
  expect(screen.getByTestId('upload-download-action')).toBeEnabled()
  expect(screen.getByTestId('upload-delete-action')).toBeDisabled()
  expect(screen.getByTestId('upload-reprocess-action')).toBeDisabled()
  expect(screen.getByTestId('upload-delete-action')).toBeDisabled()

  expect(queryVisibleForAllCheckbox()).toBeEnabled()
}

const testUnpublishedWritePermissions = async () => {
  await screen.findByTestId('logout-button')
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
  expect(screen.getByTestId('edit-metadata-button')).toBeEnabled()
  expect(screen.getByTestId('publish-upload-button')).toBeEnabled()

  expect(queryVisibleForAllCheckbox()).toBeEnabled()
}

const testEmbargoedPublishesWritePermissions = async () => {
  await screen.findByTestId('logout-button')
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if only the first two steps are shown
  expect(screen.queryByText('Prepare and upload your files')).toBeInTheDocument()
  expect(screen.queryByText('Processing completed, 1/1 entries processed')).toBeInTheDocument()
  expect(screen.queryByText('You can either select and edit individual entries from the list above, or edit all entries at once.')).toBeInTheDocument()

  expect(screen.getByTestId('edit-members-action')).toBeEnabled()
  expect(screen.getByTestId('upload-download-action')).toBeEnabled()
  expect(screen.getByTestId('upload-reprocess-action')).toBeDisabled()
  expect(screen.getByTestId('upload-delete-action')).toBeDisabled()
  expect(screen.getByTestId('edit-metadata-button')).toBeEnabled()
  expect(screen.getByText('Lift Embargo')).toBeEnabled()

  expect(queryVisibleForAllCheckbox()).toBeDisabled()
}

const testReadOnlyPermissions = async (isLoggedIn) => {
  await screen.findByTestId(isLoggedIn ? 'logout-button' : 'login-register-button')
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
  expect(screen.queryByTestId('upload-delete-action')).toBeDisabled()
  expect(screen.queryByTestId('edit-metadata-button')).not.toBeInTheDocument()
  expect(screen.queryByTestId('publish-upload-button')).not.toBeInTheDocument()

  expect(queryVisibleForAllCheckbox()).not.toBeInTheDocument()
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
  ], [
    'Published and logged in as coauthor group owner',
    'tests.states.uploads.published_coauthor_group',
    'tests/data/uploads/uploadpage-published-coauthor-group-owner',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Published and logged in as coauthor group member',
    'tests.states.uploads.published_coauthor_group',
    'tests/data/uploads/uploadpage-published-coauthor-group-member',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<><LoginLogout/><UploadPage uploadId={uploadId}/></>)
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
    'Published and logged in as reviewer group owner',
    'tests.states.uploads.published_reviewer_group',
    'tests/data/uploads/uploadpage-published-reviewer-group-owner',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Published and logged in as reviewer group member',
    'tests.states.uploads.published_reviewer_group',
    'tests/data/uploads/uploadpage-published-reviewer-group-member',
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
    'Unpublished and logged in as reviewer group owner',
    'tests.states.uploads.unpublished_reviewer_group',
    'tests/data/uploads/uploadpage-unpublished-reviewer-group-owner',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Unpublished and logged in as reviewer group member',
    'tests.states.uploads.unpublished_reviewer_group',
    'tests/data/uploads/uploadpage-unpublished-reviewer-group-member',
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
  ], [
    'Published with embargo and logged in as reviewer group owner',
    'tests.states.uploads.published_with_embargo_reviewer_group',
    'tests/data/uploads/uploadpage-published-with-embargo-reviewer-group-owner',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Published with embargo and logged in as reviewer group member',
    'tests.states.uploads.published_with_embargo_reviewer_group',
    'tests/data/uploads/uploadpage-published-with-embargo-reviewer-group-member',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<><LoginLogout/><UploadPage uploadId={uploadId}/></>)
  const isLoggedIn = username !== ''
  await testReadOnlyPermissions(isLoggedIn)
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
  ], [
    'Unpublished and logged in as coauthor group owner',
    'tests.states.uploads.unpublished_coauthor_group',
    'tests/data/uploads/uploadpage-unpublished-coauthor-group-owner',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Unpublished and logged in as coauthor group member',
    'tests.states.uploads.unpublished_coauthor_group',
    'tests/data/uploads/uploadpage-unpublished-coauthor-group-member',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<><LoginLogout/><UploadPage uploadId={uploadId}/></>)
  await testUnpublishedWritePermissions()
  await testShownColumnsAction()
})

const expectEntriesOrder = (list) => {
  const rows = screen.queryAllByTestId('datatable-row')
  const n = list.length
  expect(rows.length).toBe(n)
  for (let i = 0; i < n; i++) {
    expect(within(rows[i]).queryByText(`vasp_${list[i]}.xml`)).toBeInTheDocument()
  }
}

test('Render upload page: multiple entries', async () => {
  await startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/multiple_entries', 'test', 'password')
  render(<><LoginLogout/><UploadPage uploadId="dft_upload_1"/></>)

  await screen.findByTestId('logout-button')
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table header is rendered correctly
  expect(screen.queryByText('6 entries')).toBeInTheDocument()
  expect(screen.queryByTestId('table-pagination')).toBeInTheDocument()
  expect(screen.queryByTestId('datatable-body')).toBeInTheDocument()

  const datatableBody = screen.getByTestId('datatable-body')

  // Test the default order of the entries
  expect(within(datatableBody).queryByText('vasp_6.xml')).not.toBeInTheDocument()
  expectEntriesOrder([1, 2, 3, 4, 5])

  // Test the order of entries: ascending sort by Mainfile
  await userEvent.click(screen.getByTestId('sortable_mainfile'))
  await waitFor(() =>
    expect(within(datatableBody).queryByText('vasp_6.xml')).not.toBeInTheDocument()
  )
  expectEntriesOrder([1, 2, 3, 4, 5])

  await waitForGUI()
  // Test the order of entries: descending sort by Mainfile
  await userEvent.click(screen.getByTestId('sortable_mainfile'))
  await waitFor(() =>
    expect(within(datatableBody).queryByText('vasp_1.xml')).not.toBeInTheDocument()
  )
  expectEntriesOrder([6, 5, 4, 3, 2])
})

test('Delete selected entries from table', async () => {
  await startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/delete_entries_from_table', 'test', 'password')
  render(<><LoginLogout/><UploadPage uploadId="dft_upload_1"/></>)

  await screen.findByTestId('logout-button')
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
  let deleteConfirmButton = await screen.findByButtonText('Delete 1 entry')
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
  deleteConfirmButton = await screen.findByButtonText('Delete 2 entries')
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
  ], [
    'Published with embargo and logged in as coauthor group owner',
    'tests.states.uploads.published_with_embargo_coauthor_group',
    'tests/data/uploads/uploadpage-published-with-embargo-coauthor-group-owner',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Published with embargo and logged in as coauthor group member',
    'tests.states.uploads.published_with_embargo_coauthor_group',
    'tests/data/uploads/uploadpage-published-with-embargo-coauthor-group-member',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<><LoginLogout/><UploadPage uploadId={uploadId}/></>)
  await testEmbargoedPublishesWritePermissions()
  await testShownColumnsAction()
})

test('Toggle visible for all checkbox; check embargo, icon', async () => {
  await startAPI(
    'tests.states.uploads.unpublished',
    'tests/data/uploads/uploadpage-dialog-toggle-visible',
    'test',
    'password'
  )

  const user = userEvent.setup()

  async function testAndToggleCheckbox(initialState, { skipToggle = false } = {}) {
    const checkbox = queryVisibleForAllCheckbox()
    await waitFor(() => expect(checkbox.checked).toEqual(initialState))

    if (skipToggle) return

    user.click(checkbox)
    await waitFor(() => expect(checkbox.checked).toEqual(!initialState))
  }

  render(<><LoginLogout/><UploadPage uploadId='dft_upload'/></>)
  await screen.findByTestId('logout-button')

  // set embargo to value
  const embargoLabel = await screen.findByText('Embargo period', { selector: 'label' })
  const embargoDiv = embargoLabel.closest('div')
  const embargoButton = embargoDiv.querySelector('[role="button"]')
  const embargoHelper = within(embargoDiv).getByText('publish without embargo')
  expect(embargoButton).not.toHaveAttribute('aria-disabled', 'true')
  expect(embargoButton).toHaveTextContent('No embargo')
  await user.click(embargoButton)
  await user.click(await screen.findByRole('option', { name: '36' }))
  expect(embargoButton).toHaveTextContent('36')
  expect(embargoHelper).toHaveTextContent('months before the data becomes public')
  expect(screen.getByTooltip('Unpublished, accessible by you, coauthors and reviewers')).toBeInTheDocument()

  await testAndToggleCheckbox(false)
  expect(embargoButton).toHaveAttribute('aria-disabled', 'true')
  expect(embargoButton).toHaveTextContent('No embargo')
  expect(embargoHelper).toHaveTextContent('Upload is publicly visible, embargo disabled')
  expect(screen.getByTooltip('Unpublished but accessible by everyone')).toBeInTheDocument()

  await testAndToggleCheckbox(true)
  expect(embargoButton).not.toHaveAttribute('aria-disabled', 'true')
  expect(embargoButton).toHaveTextContent('No embargo')
  expect(embargoHelper).toHaveTextContent('publish without embargo')
  expect(screen.getByTooltip('Unpublished, accessible by you, coauthors and reviewers')).toBeInTheDocument()

  await testAndToggleCheckbox(false, { skipToggle: true })
})
