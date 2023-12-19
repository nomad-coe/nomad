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
import { renderNoAPI, screen } from '../conftest.spec'
import { fireEvent } from '@testing-library/react'
import { DateTimeEditQuantity } from './DateTimeEditQuantity'

const changeValue = async (input, newValue) => {
  await fireEvent.change(input, { target: { value: newValue } })
  await fireEvent.keyDown(input, {key: 'Enter', code: 'Enter'})
}

const handleChange = jest.fn(value => {})
const quantityDef = {
  name: 'name',
  description: `This is **MARKDOWN** help text.`
}

test('initial value is displayed correctly', async () => {
  renderNoAPI(<DateTimeEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value='1970-01-01T15:00:00.000Z'
  />)
  screen.getByDisplayValue('01/01/1970 16:00')
})

test('nil value should not display anything', async () => {
  renderNoAPI(<DateTimeEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
  />)
  const input = screen.getByRole('textbox')
  expect(input.value).toBe('')
})

test('internal change updates time correctly and triggers callback', async () => {
  renderNoAPI(<DateTimeEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
  />)
  const input = screen.getByRole('textbox')
  await changeValue(input, '01/01/1970 17:00')
  screen.getByDisplayValue('01/01/1970 17:00')
  expect(handleChange.mock.calls).toHaveLength(1)
  expect(handleChange.mock.calls[0][0]).toBe('1970-01-01T16:00:00.000Z')
})

test('external change (without internal change) updates time correctly', async () => {
  const {rerender} = renderNoAPI(<DateTimeEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
  />)
  rerender(<DateTimeEditQuantity
    quantityDef={quantityDef}
    onChange={handleChange}
    value='1970-01-01T16:00:00.000Z'
  />)
  screen.getByDisplayValue('01/01/1970 17:00')
})
