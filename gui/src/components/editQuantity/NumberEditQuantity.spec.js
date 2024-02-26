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
import {render, screen} from '../conftest.spec'
import {within} from '@testing-library/dom'
import {fireEvent, waitFor} from '@testing-library/react'
import {NumberEditQuantity} from './NumberEditQuantity'

const changeValue = (input, newValue) => {
  fireEvent.change(input, { target: { value: newValue } })
  fireEvent.keyDown(input, {key: 'Enter', code: 'Enter'})
}

const handleChange = jest.fn(value => {})

describe('Test numberEditQuantity', () => {
  it('functionality', async () => {
    render(<NumberEditQuantity
      quantityDef={{
        name: 'name',
        description: `This is **MARKDOWN** help text.`
      }}
      onChange={handleChange}
      type={{type_kind: 'python', type_data: 'int'}}
      value={10}
    />)
    const numberFieldValue = screen.queryByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('10'))

    await changeValue(numberFieldValueInput, '5')
    await changeValue(numberFieldValueInput, '')

    await waitFor(() => expect(handleChange.mock.calls).toHaveLength(2))
    await waitFor(() => expect(handleChange.mock.calls[0][0]).toBe(5))
    await waitFor(() => expect(handleChange.mock.calls[1][0]).toBe(undefined))
  })

  it('without default display unit', async () => {
    render(
      <NumberEditQuantity
        quantityDef={{
          name: 'name',
          unit: 'meter'
        }}
        type={{type_kind: 'python', type_data: 'int'}}
        value={10}
      />
    )
    const numberFieldValue = screen.queryByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    // should be rendered in default unit system
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('100000000000'))
  })

  it('with default display unit', async () => {
    render(
      <NumberEditQuantity
        quantityDef={{
          name: 'name',
          unit: 'meter',
          m_annotations: {
            display: [{
              unit: 'mm'
            }]
          }
        }}
        type={{type_kind: 'python', type_data: 'int'}}
        value={10}
      />
    )
    const numberFieldValue = screen.queryByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    // should be rendered in default unit provided by schema
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('10000'))
  })

  it('complex unit', async () => {
    render(
      <NumberEditQuantity
        quantityDef={{
          name: 'name',
          unit: 'm**2 / second**2'
        }}
        type={{type_kind: 'python', type_data: 'int'}}
        value={10}
      />
    )
    const numberFieldValue = screen.queryByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('1e-9'))
  })

  it('complex unit with display unit', async () => {
    render(
      <NumberEditQuantity
        quantityDef={{
          name: 'name',
          unit: 'm**2 / second**2',
          m_annotations: {
            display: [{
              unit: 'm**2 / second**2'
            }]
          }
        }}
        type={{type_kind: 'python', type_data: 'int'}}
        value={10}
      />
    )
    const numberFieldValue = screen.queryByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('10'))
  })

  it('with default display unit (deprecated version)', async () => {
    render(
      <NumberEditQuantity
        quantityDef={{
          name: 'name',
          unit: 'meter',
          m_annotations: {
            eln: [{
              defaultDisplayUnit: 'mm'
            }]
          }
        }}
        type={{type_kind: 'python', type_data: 'int'}}
        value={10}
      />
    )
    const numberFieldValue = screen.queryByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    // should be rendered in default unit provided by schema
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('10000'))
  })
})
