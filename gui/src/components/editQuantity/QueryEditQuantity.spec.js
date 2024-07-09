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
import {closeAPI, render, screen, startAPI, waitForGUI} from '../conftest.spec'
import QueryEditQuantity from "./QueryEditQuantity"
import {waitFor, within} from "@testing-library/dom"
import {act} from "react-dom/test-utils"
import userEvent from '@testing-library/user-event'

const handleChange = jest.fn(value => {})
const quantityDef = {
  name: 'myQuery',
  _qualifiedName: `myQuery-qualifiedName`,
  description: `
            This is **MARKDOWN** help text.
          `
}

const testSearchDialogCancelButton = async () => {
  const dialog = screen.getByTestId('search-dialog')
  await waitFor(() => expect(screen.queryByText('visibility=visible')).toBeInTheDocument())

  // cancel the search
  await userEvent.click(within(dialog).getByRole('button', {name: /cancel/i}))
  await waitFor(() => expect(screen.queryByTestId('search-dialog')).not.toBeInTheDocument())
}

const testSearchDialogOkButton = async () => {
  const dialog = screen.getByTestId('search-dialog')
  await waitFor(() => expect(screen.queryByText('visibility=visible')).toBeInTheDocument())

  // accept the search
  await userEvent.click(within(dialog).getByTestId('search-dialog-ok'))
  await waitFor(() => expect(screen.queryByTestId('search-dialog')).not.toBeInTheDocument())
}

test('Test QueryEditQuantity', async () => {
  await startAPI('tests.states.entry.eln', 'tests/data/editquantity/query', 'test', 'password')
  render(<QueryEditQuantity
    quantityDef={quantityDef}
    storeInArchive={true}
    value={{
      filters: {visibility: 'visible'},
      results: [
        {entry_id: '1', mainfile: 'a'},
        {entry_id: '2', mainfile: 'b'},
        {entry_id: '3', mainfile: 'c'}
      ]}}
    onChange={handleChange}
  />)

  const input = screen.getByRole('textbox')
  expect(input.value).toBe('')

  screen.queryByText('visibility:visible')
  screen.queryByText('3 results')

  const searchDialogButton = screen.getByTitle('Search dialog').closest('button')
  expect(searchDialogButton).toBeEnabled()

  await act(async () => { userEvent.click(searchDialogButton) })
  await waitForGUI(1000, true)

  await testSearchDialogCancelButton()
  screen.getByText('visibility:visible')
  expect(input).toHaveAttribute('placeholder', '3 results')

  await act(async () => { userEvent.click(searchDialogButton) })
  await waitForGUI(1000, true)

  await testSearchDialogOkButton()
  screen.getByText('visibility:visible')

  // Assert the new results
  await waitFor(() => expect(handleChange.mock.calls[0][0].data[0].entry_id).toBe('bC7byHvWJp62Sn9uiuJUB38MT5j-'))
  await waitFor(() => expect(handleChange.mock.calls[0][0].data[0].mainfile).toBe('sample.archive.json'))
  await waitFor(() => expect(handleChange.mock.calls[0][0].data[1].entry_id).toBe('83DS7AzwqTKFVwlrdVeaL3kMSLU_'))
  await waitFor(() => expect(handleChange.mock.calls[0][0].data[1].mainfile).toBe('schema.archive.yaml'))

  // Clear results
  const clearResultsButton = screen.getByTitle('Clear results').closest('button')
  expect(clearResultsButton).toBeEnabled()
  await act(async () => { userEvent.click(clearResultsButton) })

  await waitFor(() => expect(handleChange.mock.calls[1][0]).toBe(undefined))

  closeAPI()
})
