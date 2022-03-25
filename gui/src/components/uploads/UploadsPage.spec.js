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
import {UploadsPage} from './UploadsPage'
import {withLoginRequired} from '../api'
import {within} from '@testing-library/dom'
import {fireEvent, waitFor} from '@testing-library/react'

test('Render uploads page: not authenticated', async () => {
  startAPI('tests.states.uploads.published', 'tests/data/uploads/uploads_not_authenticated')
  let Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('You have to login to use this functionality.')
  await expect(screen.getByTestId('login-register-button')).toBeEnabled()

  closeAPI()
})

test('Render uploads page: authenticated', async () => {
  startAPI('tests.states.uploads.empty', 'tests/data/uploads/empty', 'test', 'password')
  let Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('Create a new upload')
  expect(screen.queryByText('Your existing uploads')).not.toBeInTheDocument()

  closeAPI()
})

test('Render uploads page: maximum unpublished exceeded', async () => {
  startAPI('tests.states.uploads.maximum_unpublished', 'tests/data/uploads/maximum_unpublished', 'test', 'password')
  let Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('Create a new upload')
  expect(screen.queryByText('Your existing uploads')).toBeInTheDocument()
  expect(screen.queryByText('Create a new upload')).toBeDisabled()
  expect(screen.queryByText('You have reached maximum number of unpublished uploads!')).toBeInTheDocument()

  closeAPI()
})

test('Render uploads page: sort by upload create time', async () => {
  startAPI('tests.states.uploads.multiple_uploads', 'tests/data/uploads/multiple_uploads', 'test', 'password')
  let Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await waitFor(() => {
    expect(screen.queryByText('1-10 of 11')).toBeInTheDocument()
  })

  expect(screen.queryByText('Your existing uploads')).toBeInTheDocument()
  expect(screen.queryByText('Create a new upload')).toBeEnabled()
  expect(screen.queryByText('You have reached maximum number of unpublished uploads!')).not.toBeInTheDocument()

  expect(screen.queryByTestId('table-pagination')).toBeInTheDocument()
  expect(screen.queryByTestId('datatable-body')).toBeInTheDocument()

  let datatableBody = screen.getByTestId('datatable-body')

  // Test if the pagination works correctly
  let rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(10)
  expect(within(datatableBody).queryByText('dft_upload_1')).not.toBeInTheDocument()

  // Test the order of uploads: by default is descending upload create time
  for (let i = 0; i < 10; i++) {
    expect(within(rows[i]).queryByText(`dft_upload_${11 - i}`)).toBeInTheDocument()
    expect(within(rows[i]).queryByTitle(((i + 1) % 2 === 0 ? 'published upload' : 'this upload is not yet published'))).toBeInTheDocument()
  }

  // Test the order of uploads: ascending sort by upload create time
  fireEvent.click(screen.queryByTestId('sortable_upload_create_time'))
  rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(10)
  await waitFor(() =>
    expect(within(rows[9]).queryByText(`dft_upload_10`)).toBeInTheDocument()
  )

  expect(within(datatableBody).queryByText('dft_upload_11')).not.toBeInTheDocument()
  for (let i = 0; i < 10; i++) {
    expect(within(rows[i]).queryByText(`dft_upload_${i + 1}`)).toBeInTheDocument()
    expect(within(rows[i]).queryByTitle(((i + 1) % 2 === 0 ? 'published upload' : 'this upload is not yet published'))).toBeInTheDocument()
  }

  closeAPI()
})

test('Render uploads page: sort by upload name', async () => {
  startAPI('tests.states.uploads.multiple_uploads', 'tests/data/uploads/sort_by_upload_name', 'test', 'password')
  let Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await waitFor(() => {
    expect(screen.queryByText('1-10 of 11')).toBeInTheDocument()
  })

  // Test the order of uploads: sort by Upload name
  fireEvent.click(screen.queryByTestId('sortable_upload_name'))
  let rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(10)
  await waitFor(() =>
    expect(within(rows[9]).queryByText(`dft_upload_10`)).toBeInTheDocument()
  )
  for (let i = 0; i < 10; i++) expect(within(rows[i]).queryByText(`dft_upload_${i + 1}`)).toBeInTheDocument()

  closeAPI()
})

test('Render uploads page: sort by published', async () => {
  startAPI('tests.states.uploads.multiple_uploads', 'tests/data/uploads/sort_by_published', 'test', 'password')
  let Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await waitFor(() => {
    expect(screen.queryByText('1-10 of 11')).toBeInTheDocument()
  })

  // Test the order of uploads: sort by Published
  fireEvent.click(screen.queryByTestId(`sortable_published`))
  let rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(10)
  await waitFor(() =>
    expect(within(rows[9]).queryByText(`dft_upload_8`)).toBeInTheDocument()
  )
  for (let i = 0; i < 6; i++) expect(within(rows[i]).queryByTitle('this upload is not yet published')).toBeInTheDocument()
  for (let i = 6; i < 10; i++) expect(within(rows[i]).queryByTitle('published upload')).toBeInTheDocument()
  expect(within(rows[0]).queryByText(`dft_upload_1`)).toBeInTheDocument()
  expect(within(rows[1]).queryByText(`dft_upload_3`)).toBeInTheDocument()
  expect(within(rows[2]).queryByText(`dft_upload_5`)).toBeInTheDocument()
  expect(within(rows[3]).queryByText(`dft_upload_7`)).toBeInTheDocument()
  expect(within(rows[4]).queryByText(`dft_upload_9`)).toBeInTheDocument()
  expect(within(rows[5]).queryByText(`dft_upload_11`)).toBeInTheDocument()
  expect(within(rows[6]).queryByText(`dft_upload_2`)).toBeInTheDocument()
  expect(within(rows[7]).queryByText(`dft_upload_4`)).toBeInTheDocument()
  expect(within(rows[8]).queryByText(`dft_upload_6`)).toBeInTheDocument()

  closeAPI()
})
