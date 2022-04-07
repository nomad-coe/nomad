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
import { fireEvent } from '@testing-library/dom'
import { startAPI, closeAPI, screen } from '../../conftest.spec'
import { renderSearchEntry, expectInputHeader } from '../conftest.spec'
import userEvent from '@testing-library/user-event'
import InputSlider from './InputSlider'

const min = '1'
const max = '5'
const stateName = 'tests.states.search.search'

describe('', () => {
  const quantity = 'results.material.n_elements'
  const step = 1
  beforeEach(() => {
    // API state with single min_max aggregation result
    startAPI(stateName, 'tests/data/search/min_max_aggregation_n_elements')

    // Render InputSlider within an entry search context.
    renderSearchEntry(<InputSlider quantity={quantity} step={step}/>)
  })
  afterEach(() => closeAPI())

  test('initial state is loaded correctly', async () => {
    // Test immediately displayed elements
    expectInputHeader(quantity, true)
    expect(screen.getByText('min')).toBeInTheDocument()
    expect(screen.getByText('max')).toBeInTheDocument()

    // Test elements that are displayed after API response
    expect(await screen.findByDisplayValue(min)).toBeInTheDocument()
    expect(await screen.findByDisplayValue(max)).toBeInTheDocument()
  })

  describe('an error message is displayed when', () => {
    test.each([
      ['minimum', 'non-numeric'],
      ['minimum', 'empty'],
      ['maximum', 'non-numeric'],
      ['maximum', 'empty']
    ])(
      '%s field is %s',
      async (fieldType, inputType) => {
        const displayValue = fieldType === 'maximum' ? max : min
        const input = {
          'non-numeric': 'hello',
          'empty': ' '
        }[inputType]
        const field = await screen.findByDisplayValue(displayValue)
        userEvent.type(field, input)
        fireEvent.keyDown(field, {key: 'Enter', code: 'Enter'})
        expect(screen.getByText(`Invalid ${fieldType} value.`)).toBeInTheDocument()
      }
    )
  })

  describe('the following numeric values are accepted:', () => {
    test.each([
      ['1'],
      ['1.0'],
      ['-1.0'],
      ['-1.0e5'],
      ['-1.0e-5']
    ])(
      '%s',
      async (input) => {
        const field = await screen.findByDisplayValue(min)
        userEvent.type(field, input)
        fireEvent.keyDown(field, {key: 'Enter', code: 'Enter'})
        expect(screen.queryByText(`Invalid minimum value.`)).toBeNull()
      }
    )
  })

  test('slider positions are correctly set', async () => {
    const inputMin = await screen.findByDisplayValue(min)
    const inputMax = await screen.findByDisplayValue(max)
    const sliders = screen.getAllByRole('slider')
    const sliderMin = sliders[0]
    const sliderMax = sliders[1]

    // Initially both sliders at the ends
    expect(sliderMin).toHaveStyle(`left: 0%`)
    expect(sliderMax).toHaveStyle(`left: 100%`)

    // After changing the min field, min slider moves
    userEvent.type(inputMin, '3')
    fireEvent.keyDown(inputMin, {key: 'Enter', code: 'Enter'})
    expect(sliderMin).toHaveStyle(`left: 50%`)

    // After changing the max field, max slider moves
    userEvent.type(inputMax, '3')
    fireEvent.keyDown(inputMax, {key: 'Enter', code: 'Enter'})
    expect(sliderMax).toHaveStyle(`left: 50%`)
  })

  test('changing sliders changes the input values', async () => {
    const inputMin = await screen.findByDisplayValue(min)
    const inputMax = await screen.findByDisplayValue(max)
    const sliders = screen.getAllByRole('slider')
    const sliderMin = sliders[0]
    const sliderMax = sliders[1]

    // Moving min slider changes min field
    fireEvent.keyDown(sliderMin, {key: 'Up', code: 'Up'})
    expect(inputMin.value).toBe((+min + step).toString())

    // Moving max slider changes max field
    fireEvent.keyDown(sliderMax, {key: 'Down', code: 'Down'})
    expect(inputMax.value).toBe((+max - step).toString())
  })
})
