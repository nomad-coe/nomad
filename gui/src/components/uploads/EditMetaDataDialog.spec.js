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
import {fireEvent, waitFor} from '@testing-library/react'

const testComment = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Open the members dialog
  fireEvent.click(screen.getByTestId('edit-metadata-button'))
  await waitFor(() => expect(screen.queryByText('Submit')).toBeInTheDocument())

  expect(screen.queryByButtonText('Submit')).toBeDisabled()
  const dialog = screen.getByTestId('edit-metadata-dialog')
  expect(within(dialog).queryByText('Edit upload meta data')).toBeInTheDocument()
  expect(within(dialog).queryByText('Comments')).toBeInTheDocument()
  expect(within(dialog).queryByText('Mocked')).toBeInTheDocument()

  const commentField = within(dialog).getByTestId('metadata-comment-field')
  fireEvent.change(commentField, {target: {value: 'new comment'}})
  await waitFor(() => expect(screen.queryByText('Submit')).toBeEnabled())
}

const testReferences = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Open the members dialog
  fireEvent.click(screen.getByTestId('edit-metadata-button'))
  await waitFor(() => expect(screen.queryByText('Submit')).toBeInTheDocument())

  expect(screen.queryByButtonText('Submit')).toBeDisabled()
  const dialog = screen.getByTestId('edit-metadata-dialog')
  expect(within(dialog).queryByText('Edit upload meta data')).toBeInTheDocument()
  expect(within(dialog).queryByText('References')).toBeInTheDocument()
  expect(within(dialog).queryByText('Datasets')).toBeInTheDocument()
  expect(within(dialog).queryByTestId('reference-add-button')).toBeDisabled()

  let rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(1)

  expect(within(rows[0]).queryByText('doi')).toBeInTheDocument()
  expect(within(rows[0]).queryByTestId('reference-delete-action')).toBeEnabled()

  fireEvent.click(within(rows[0]).getByTestId('reference-delete-action'))
  await waitFor(() => expect(screen.queryByText('doi')).not.toBeInTheDocument())
  rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(0)
  expect(screen.queryByTestId('reference-delete-action')).not.toBeInTheDocument()

  const referenceField = within(dialog).getByTestId('new-reference-field')
  const referenceAddButton = within(dialog).queryByTestId('reference-add-button')

  fireEvent.change(referenceField, {target: {value: 'invalid url'}})
  await waitFor(() => expect(referenceAddButton).toBeEnabled())

  fireEvent.click(referenceAddButton)
  await waitFor(() => expect(within(dialog).queryByText('Pleas enter a valid URL ...')).toBeInTheDocument())
  rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(0)

  fireEvent.change(referenceField, {target: {value: 'http://valid.url'}})
  await waitFor(() => expect(referenceAddButton).toBeEnabled())

  fireEvent.click(referenceAddButton)
  await waitFor(() => expect(within(dialog).queryByText('Pleas enter a valid URL ...')).not.toBeInTheDocument())
  rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(1)

  expect(within(rows[0]).queryByText('http://valid.url')).toBeInTheDocument()
  expect(within(rows[0]).queryByTestId('reference-delete-action')).toBeEnabled()
  expect(screen.queryByText('Submit')).toBeEnabled()
}

const testDatasetForPublished = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Open the members dialog
  fireEvent.click(screen.getByTestId('edit-metadata-button'))
  await waitFor(() => expect(screen.queryByText('Submit')).toBeInTheDocument())

  expect(screen.queryByButtonText('Submit')).toBeDisabled()
  const dialog = screen.getByTestId('edit-metadata-dialog')
  expect(within(dialog).queryByText('Edit upload meta data')).toBeInTheDocument()
  expect(within(dialog).queryByText('Comments')).toBeInTheDocument()
  expect(within(dialog).queryByText('Mocked')).toBeInTheDocument()

  fireEvent.click(screen.getByTestId('reference-delete-action'))
  await waitFor(() => expect(screen.queryByText('doi')).not.toBeInTheDocument())

  const datasetField = within(dialog).getByTestId('new-dataset-field')
  fireEvent.change(datasetField, {target: {value: 'new dataset (1)'}})
  await waitFor(() => expect(within(dialog).queryByText('add entry to new dataset')).toBeEnabled())

  fireEvent.click(within(dialog).queryByText('add entry to new dataset'))
  await waitFor(() => expect(within(dialog).queryByText('add entry to new dataset')).not.toBeInTheDocument())
  const rows = within(dialog).queryAllByTestId('datatable-row')
  expect(rows.length).toBe(1)
  await waitFor(() => expect(within(dialog).queryByText('Open in new tab')).not.toBeInTheDocument())

  expect(within(rows[0]).queryByText('new dataset (1)')).toBeInTheDocument()
  expect(within(rows[0]).queryByTestId('dataset-delete-action')).toBeEnabled()
  expect(screen.queryByText('Submit')).toBeEnabled()
}

const testReadOnlyPermissions = async () => {
  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Members dialog should be disabled
  expect(screen.queryByTestId('edit-metadata-button')).not.toBeInTheDocument()
}

test.each([
  [
    'Published and logged in as main author',
    'tests.states.uploads.published',
    'tests/data/uploads/metadata-dialog-published-author',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published and logged in as coauthor',
    'tests.states.uploads.published',
    'tests/data/uploads/metadata-dialog-published-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Unpublished and logged in as main author',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/metadata-dialog-unpublished-author',
    'dft_upload',
    'test',
    'password'
  ], [
    'Unpublished and logged in as coauthor',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/metadata-dialog-unpublished-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ], [
    'Published with embargo and logged in as main author',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/metadata-dialog-published-with-embargo-author',
    'dft_upload',
    'test',
    'password'
  ], [
    'Published with embargo and logged in as coauthor',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/metadata-dialog-published-with-embargo-coauthor',
    'dft_upload',
    'scooper',
    'password'
  ]
])('Metadata dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testComment()
  await testReferences()
  await testDatasetForPublished()
  closeAPI()
})

test.each([
  [
    'Published and logged in as reviewer',
    'tests.states.uploads.published',
    'tests/data/uploads/metadata-dialog-published-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published and logged in as neither reviewer nor coauthor or main author',
    'tests.states.uploads.published',
    'tests/data/uploads/metadata-dialog-published-external',
    'dft_upload',
    'admin',
    'password'
  ], [
    'Published and not authenticated',
    'tests.states.uploads.published',
    'tests/data/uploads/metadata-dialog-published-noauth',
    'dft_upload',
    '',
    ''
  ], [
    'Unpublished and logged in as reviewer',
    'tests.states.uploads.unpublished',
    'tests/data/uploads/metadata-dialog-unpublished-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ], [
    'Published with embargo and logged in as reviewer',
    'tests.states.uploads.published_with_embargo',
    'tests/data/uploads/metadata-dialog-published-with-embargo-reviewer',
    'dft_upload',
    'ttester',
    'password'
  ]
])('Metadata dialog: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)
  await testReadOnlyPermissions()
  closeAPI()
})
