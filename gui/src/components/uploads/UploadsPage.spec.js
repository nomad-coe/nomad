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
import {exampleUploads as exampleUploadsJSON} from '../../config'
import {UploadsPage} from './UploadsPage'
import {withLoginRequired} from '../api'
import {within} from '@testing-library/dom'
import {fireEvent, waitFor, waitForElementToBeRemoved} from '@testing-library/react'

test('Render uploads page: sort by upload create time', async () => {
  await startAPI('tests.states.uploads.multiple_uploads', 'tests/data/uploads/uploadspage', 'test', 'password')
  const Component = withLoginRequired(UploadsPage, undefined)
  render(<Component/>)

  // Wait to load the page, i.e. wait for some text to appear
  await waitFor(() => {
    expect(screen.queryByText('1-10 of 11')).toBeInTheDocument()
  })

  const datatableBody = screen.getByTestId('datatable-body')

  // Test if the pagination works correctly
  let rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(10)
  expect(within(datatableBody).queryByText('dft_upload_1')).not.toBeInTheDocument()

  // Test the order of uploads: by default is descending upload create time
  for (let i = 0; i < 10; i++) {
    expect(within(rows[i]).queryByText(`dft_upload_${11 - i}`)).toBeInTheDocument()
    expect(within(rows[i]).queryByTitle(((i + 1) % 2 === 0
      ? 'Published and accessible by everyone'
      : 'Unpublished, only accessible by you, coauthors and reviewers'
    ))).toBeInTheDocument()
  }

  // Test the order of uploads: sort by Upload name
  fireEvent.click(screen.queryByTestId('sortable_upload_name'))
  rows = screen.queryAllByTestId('datatable-row')
  expect(rows.length).toBe(10)
  await waitFor(() =>
    expect(within(rows[9]).queryByText(`dft_upload_10`)).toBeInTheDocument()
  )
  for (let i = 0; i < 10; i++) expect(within(rows[i]).queryByText(`dft_upload_${i + 1}`)).toBeInTheDocument()

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
    expect(within(rows[i]).queryByTitle(((i + 1) % 2 === 0
      ? 'Published and accessible by everyone'
      : 'Unpublished, only accessible by you, coauthors and reviewers'
    ))).toBeInTheDocument()
  }

  // Testing the "Add example uploads" functionality
  const exampleButton = screen.getByRole('button', { name: /add example uploads/i })

  // Testing for Example Button to be in the document
  expect(exampleButton).toBeInTheDocument()

  // Testing for the example button to show the Dialog box with a list of examples
  fireEvent.click(exampleButton)
  const dialogTitle = screen.getByRole('heading', { name: /select a sample upload/i })
  expect(dialogTitle).toBeInTheDocument()

  // Testing the example list to contain the exampleUploads json file
  const exampleUploads = Object.keys(exampleUploadsJSON).reduce((uploads, category) => {
    const categoryJSON = exampleUploadsJSON[category]
    return [...uploads, ...Object.keys(categoryJSON).map(upload => categoryJSON[upload])]
  }, [])
  exampleUploads.forEach(upload => {
    expect(screen.queryByText(upload.title)).toBeInTheDocument()
  })

  // Testing for all add-buttons to be present for each example
  expect(screen.queryAllByRole('button', { name: /add/i }).length).toBe(exampleUploads.length)

  // Testing for the cancel button in the dialog to close it
  fireEvent.click(screen.getByRole('button', { name: /cancel/i }))
  waitForElementToBeRemoved(dialogTitle)

  closeAPI()
})
