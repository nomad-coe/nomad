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
import {EditQuantityExamples} from './EditQuantityExamples'
import {within} from '@testing-library/dom'
import {fireEvent, waitFor} from '@testing-library/react'

test('correctly renders edit quantities', async () => {
  closeAPI()
  await startAPI('tests.states.entry.dft', 'tests/data/entry/dft')
  render(<EditQuantityExamples />)

  // Wait to load the entry metadata, i.e. wait for some texts to appear
  await screen.findByText('three')

  const numberFieldValue = screen.queryAllByTestId('number-edit-quantity-value')
  const numberFieldValueInput = within(numberFieldValue[0]).getByRole('textbox')
  fireEvent.change(numberFieldValueInput, { target: { value: '1.5' } })
  await waitFor(() => expect(numberFieldValueInput.value).toEqual('1.5'))

  const numberFieldUnit = screen.queryAllByTestId('number-edit-quantity-unit')
  const numberFieldUnitButton = within(numberFieldUnit[0]).getByRole('button')
  fireEvent.click(numberFieldUnitButton)
  // fireEvent.change(numberFieldUnit[0], { target: { value: 'm' } })
  // await waitFor(() => expect(numberFieldValueInput.value).toEqual('7.937658163559664e-11'))

  closeAPI()
})
