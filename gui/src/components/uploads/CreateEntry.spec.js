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

const testSectionSelectAutocomplete = async () => {
  await waitFor(() => expect(screen.queryAllByTestId('section-select-entry-activated').length).toBe(5))

  const sectionSelectEntries = screen.getAllByTestId('section-select-entry-activated')

  await waitFor(() => expect(within(sectionSelectEntries[0]).queryByText('Set your name here')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[0]).queryByText('upload id: references_upload_id1')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[1]).queryByText('Set your name here')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[1]).queryByText('upload id: references_upload_id1')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[2]).queryByText('Set your name here')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[2]).queryByText('upload id: references_upload_id2')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[3]).queryByText('Set your name here')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[3]).queryByText('upload id: references_upload_id2')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[4]).queryByText('Set your name here')).toBeInTheDocument())
  await waitFor(() => expect(within(sectionSelectEntries[4]).queryByText('upload id: references_upload_id2')).toBeInTheDocument())
}

test.each([
  [
    'Test external upload referencing',
    'tests.states.entry.references',
    'tests/data/uploads/external-upload-references',
    'references_upload_id2',
    'test',
    'password'
  ]
])('Upload page: %s', async (name, state, snapshot, uploadId, username, password) => {
  await startAPI(state, snapshot, username, password)
  render(<UploadPage uploadId={uploadId}/>)

  await waitFor(() => expect(screen.getByText('Processing completed, 3/3 entries processed')).toBeInTheDocument())

  const createEntryButton = screen.getByButtonText('Create from schema')
  await act(async () => { await userEvent.click(createEntryButton) })

  await waitFor(() => expect(screen.queryByText('Select a schema')).toBeInTheDocument())
  const dialog = screen.getByTestId('create-entry-dialog')

  expect(within(dialog).queryByTestId('custom-select-schema')).not.toBeInTheDocument()
  const builtinField = within(dialog).getByTestId('builtin-select-schema')
  const builtinFieldInput = within(builtinField).getByRole('textbox')

  await waitFor(() => expect(builtinFieldInput.value).not.toEqual(''))

  const customSchemaRadio = within(dialog).getByTestId('custom-schema-radio')
  await act(async () => { await userEvent.click(customSchemaRadio) })

  expect(within(dialog).queryByTestId('builtin-select-schema')).not.toBeInTheDocument()
  const customField = within(dialog).getByTestId('custom-select-schema')
  const customFieldInput = within(customField).getByRole('textbox')

  await waitFor(() => expect(customFieldInput.value).toEqual(''))

  act(() => { fireEvent.change(customFieldInput, { target: { value: 'set' } }) })
  await waitForGUI(2000) // the input changes have been debounced

  await testSectionSelectAutocomplete()

  closeAPI()
})
