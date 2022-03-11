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
import {
  render,
  screen,
  startAPI,
  closeAPI
} from '../conftest'
import UploadPage from './UploadPage'
import { within } from '@testing-library/dom'

test('Render upload page: published|not reader|not writer', async () => {
  startAPI('tests.states.uploads.published', 'tests/data/uploads/published')
  render(<UploadPage uploadId={'dft_upload'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  expect(screen.getByTestId('edit-members-action')).toBeDisabled()
  expect(screen.getByTestId('upload-download-action')).toBeEnabled()
  expect(screen.getByTestId('upload-reprocess-action')).toBeDisabled()
  expect(screen.getByTestId('source-api-action')).toBeEnabled()
  expect(screen.getByTestId('upload-delete-action')).toBeDisabled()

  closeAPI()
})

test('Render upload page: unpublished|not reader|not writer', async () => {
  startAPI('tests.states.uploads.unpublished', 'tests/data/uploads/unpublished')
  render(<UploadPage uploadId={'dft_upload'}/>)

  await screen.findByText('You do not have access to the specified upload - not published yet.')

  closeAPI()
})

test('Render upload page: published under 3 months embargo/not reader/not writer', async () => {
  startAPI('tests.states.uploads.published_with_embargo', 'tests/data/uploads/published_with_embargo')
  render(<UploadPage uploadId={'dft_upload'}/>)

  await screen.findByText('You do not have access to the specified upload - published with embargo.')

  closeAPI()
})

test('Render upload page: multiple entries', async () => {
  startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/multiple_entries')
  render(<UploadPage uploadId={'dft_upload_1'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  expect(screen.queryByTestId('step-prepare-and-upload-your-files')).toBeInTheDocument()
  expect(screen.queryByTestId('step-process-data')).toBeInTheDocument()
  expect(screen.queryByTestId('step-edit-author-metadata')).toBeNull()
  expect(screen.queryByTestId('step-publish')).toBeNull()

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

test('Render upload page: one entry', async () => {
  startAPI('tests.states.uploads.published', 'tests/data/uploads/published')
  render(<UploadPage uploadId={'dft_upload'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table title is rendered correctly
  expect(screen.queryByText('1 entry')).toBeInTheDocument()

  closeAPI()
})

test('Render upload page: not exists', async () => {
  startAPI('tests.states.uploads.published', 'tests/data/uploads/not_exists')
  render(<UploadPage uploadId={'a_not_exists_upload_ID'}/>)

  await screen.findByText('The specified upload_id was not found.')

  closeAPI()
})
