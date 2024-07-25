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
import {renderNoAPI, screen} from '../conftest.spec'
import {BoolEditQuantity} from './BoolEditQuantity'
import {fireEvent} from "@testing-library/react"
import {act} from 'react-dom/test-utils'
import userEvent from '@testing-library/user-event'

const handleChange = jest.fn()
const quantityDef = {
  name: 'name',
  description: `This is **MARKDOWN** help text.`
}
const booleanLabels = {true: 'On', false: 'Off', undefined: 'Unknown'}

const getRadioElements = () => {
  const radioOn = screen.getByLabelText(booleanLabels[true])
  const radioOff = screen.getByLabelText(booleanLabels[false])
  const radioUnknown = screen.getByLabelText(booleanLabels[undefined])
  return {radioOn, radioOff, radioUnknown}
}

test('checkbox initial value is displayed correctly', async () => {
  renderNoAPI(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={true}
  />)
  const checkbox = screen.getByRole('checkbox')
  expect(checkbox.checked).toBe(true)
  expect(checkbox.indeterminate).toBe(false)
})

test('radio buttons are displayed correctly with initial value', async () => {
  renderNoAPI(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={true}
    booleanLabels={booleanLabels}
  />)
  const {radioOn, radioOff, radioUnknown} = getRadioElements()
  expect(radioOn).toBeInTheDocument()
  expect(radioOff).toBeInTheDocument()
  expect(radioUnknown).toBeInTheDocument()
  expect(radioOn.checked).toBe(true)
  expect(radioOff.checked).toBe(false)
  expect(radioUnknown.checked).toBe(false)
})

test('external change (without internal change) updates correctly for checkbox', async () => {
  const {rerender} = renderNoAPI(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={false}
  />)
  rerender(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={true}
  />)
  const checkbox = screen.getByRole('checkbox')
  expect(checkbox.checked).toBe(true)
  expect(checkbox.indeterminate).toBe(false)
})

test('external change (without internal change) updates correctly for radio buttons', async () => {
  const {rerender} = renderNoAPI(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={false}
    booleanLabels={booleanLabels}
  />)
  rerender(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={true}
    booleanLabels={booleanLabels}
  />)
  const {radioOn, radioOff, radioUnknown} = getRadioElements()
  expect(radioOn.checked).toBe(true)
  expect(radioOff.checked).toBe(false)
  expect(radioUnknown.checked).toBe(false)
})

test('handleChange gets called with correct arguments when user clicks on checkbox', async () => {
  renderNoAPI(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={false}
  />)
  const checkbox = screen.getByRole('checkbox')
  fireEvent.click(checkbox)
  expect(handleChange).toHaveBeenCalledWith(true)
})

test('handleChange gets called with correct arguments when user clicks on radio buttons', async () => {
  renderNoAPI(<BoolEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value={false}
    booleanLabels={booleanLabels}
  />)
  const {radioOn, radioUnknown} = getRadioElements()
  await act(async () => {
    await userEvent.click(radioOn)
  })
  expect(handleChange).toHaveBeenCalledWith(true)
  await act(async () => {
    await userEvent.click(radioUnknown)
  })
  expect(handleChange).toHaveBeenCalledWith(undefined)
})
