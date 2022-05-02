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
  screen, wait
} from '../conftest.spec'
import {EditQuantityExamples} from './EditQuantityExamples'
import {within} from '@testing-library/dom'
import {fireEvent, waitFor} from '@testing-library/react'

test('correctly renders edit quantities', async () => {
  render(<EditQuantityExamples />)

  // Wait to load the entry metadata, i.e. wait for some texts to appear
  await screen.findByText('three')

  const numberFieldValue = screen.queryAllByTestId('number-edit-quantity-value')
  const numberFieldValueInput = within(numberFieldValue[1]).getByRole('textbox')
  const numberFieldUnit = screen.queryAllByTestId('number-edit-quantity-unit')
  const numberFieldUnitInput = within(numberFieldUnit[0]).getByRole('textbox', { hidden: true })
  const numberFieldCheckBoxs = screen.queryAllByTitle('If checked, numeric value is converted when the unit is changed.')
  const numberFieldCheckBox = within(numberFieldCheckBoxs[0]).getByRole('checkbox', { hidden: true })
  await waitFor(() => expect(numberFieldUnitInput.value).toEqual('bohr'))

  // First test read mode
  await waitFor(() => expect(numberFieldCheckBox.checked).toBe(true))

  fireEvent.change(numberFieldValueInput, { target: { value: '1.5' } })
  await waitFor(() => expect(screen.queryByText(/"float_unit": 7\.937658163559664e-11/i)).toBeInTheDocument())

  fireEvent.change(numberFieldUnitInput, { target: { value: 'angstrom' } })
  await waitFor(() => expect(numberFieldValueInput.value).toEqual('0.7937658163559663'))
  await waitFor(() => expect(screen.queryByText(/"float_unit": 7\.937658163559664e-11/i)).toBeInTheDocument())

  fireEvent.change(numberFieldValueInput, { target: { value: '1.5m' } })
  await waitFor(() => expect(numberFieldValueInput.value).toEqual('1.5'))
  await waitFor(() => expect(numberFieldUnitInput.value).toEqual('meter'))
  await waitFor(() => expect(screen.queryByText(/"float_unit": 1\.5/i)).toBeInTheDocument())

  // now test write mode
  fireEvent.click(numberFieldCheckBox)
  await wait(undefined, 100)
  await waitFor(() => expect(numberFieldCheckBox.checked).toBe(false))

  fireEvent.change(numberFieldUnitInput, { target: { value: 'bohr' } })
  await waitFor(() => expect(screen.queryByText(/"float_unit": 7\.937658163559664e-11/i)).toBeInTheDocument())
  await waitFor(() => expect(numberFieldValueInput.value).toEqual('1.5'))

  fireEvent.change(numberFieldValueInput, { target: { value: '1.5angstrom' } })
  await waitFor(() => expect(numberFieldValueInput.value).toEqual('1.5'))
  await waitFor(() => expect(numberFieldUnitInput.value).toEqual('angstrom'))
  await waitFor(() => expect(screen.queryByText(/"float_unit": 1\.5e-10/i)).toBeInTheDocument())
})
