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

test('Render upload page: unauthenticated', async () => {
  startAPI('tests.states.uploads.multiple_entries', 'tests/data/uploads/multiple_entries')
  render(<UploadPage uploadId={'dft_upload_1'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  expect(screen.queryByTestId('step-prepare-and-upload-your-files')).toBeInTheDocument()
  expect(screen.queryByTestId('step-process-data')).toBeInTheDocument()
  expect(screen.queryByTestId('step-edit-author-metadata')).toBeNull()
  expect(screen.queryByTestId('step-publish')).toBeNull()

  // Test if the table header is rendered correctly
  expect(screen.queryByText('4 entries')).toBeInTheDocument()
  expect(screen.queryByRole('table-pagination')).toBeInTheDocument()
  expect(screen.queryByRole('datatable-body')).toBeInTheDocument()

  let datatableBody = screen.getByRole('datatable-body')

  // Test if the name of the entries are rendered
  expect(within(datatableBody).queryByText('vasp_1.xml')).toBeInTheDocument()
  expect(within(datatableBody).queryByText('vasp_2.xml')).toBeInTheDocument()
  expect(within(datatableBody).queryByText('vasp_3.xml')).toBeInTheDocument()
  expect(within(datatableBody).queryByText('vasp_4.xml')).toBeInTheDocument()

  closeAPI()
})

test('Render upload page: one entry', async () => {
  startAPI('tests.states.uploads.one_entry', 'tests/data/uploads/one_entry')
  render(<UploadPage uploadId={'dft_upload'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  // Test if the table header is rendered correctly
  expect(screen.queryByText('1 entry')).toBeInTheDocument()

  closeAPI()
})

test('Render upload page: not exists', async () => {
  startAPI('tests.states.uploads.one_entry', 'tests/data/uploads/not_exists')
  render(<UploadPage uploadId={'a_not_exists_upload_ID'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('The specified upload_id was not found.')

  closeAPI()
})
