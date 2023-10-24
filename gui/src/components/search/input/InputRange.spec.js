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
import userEvent from '@testing-library/user-event'
import { format } from 'date-fns'
import { startAPI, closeAPI, screen } from '../../conftest.spec'
import { renderSearchEntry, expectInputRange } from '../conftest.spec'
import InputRange from './InputRange'
import { defaultFilterData } from '../FilterRegistry'
import { DType, formatNumber } from '../../../utils'

const nBins = 30
const stateName = 'tests.states.search.histograms'
const discrete = ['results.material.n_elements', false, 1, 5, 1]
const discrete_histogram = ['results.material.n_elements', true, 1, 5, 1]
const continuous = ['results.properties.electronic.band_structure_electronic.band_gap.value', false, 0, 2, 0.1]
const continuous_histogram = ['results.properties.electronic.band_structure_electronic.band_gap.value', true, 0, 2, 2 / nBins]
const time = ['upload_create_time', false, 1585872000000, 1585872240000, (1585872240000 - 1585872000000) / nBins]
const time_histogram = ['upload_create_time', true, 1585872000000, 1585872240000, undefined]

describe('test initial state', () => {
  beforeAll(async () => { await startAPI(stateName, 'tests/data/search/inputrange-init') })
  afterAll(() => closeAPI())

  test.each([
    continuous,
    continuous_histogram,
    discrete,
    discrete_histogram,
    time,
    time_histogram
  ])('quantity: %s, histogram: %s', async (quantity, histogram, min, max) => {
    renderSearchEntry(<InputRange visible quantity={quantity} disableHistogram={!histogram}/>)
    await expectInputRange(quantity, false, histogram, false, min, max)
  })
})

describe('test invalid/valid numeric input', () => {
  beforeAll(async () => { await startAPI(stateName, 'tests/data/search/inputrange-validation') })
  afterAll(() => closeAPI())
  const [quantity, histogram, min, max] = discrete

  for (const isMin of [false, true]) {
    const value = isMin ? min : max
    const field = isMin ? 'minimum' : 'maximum'
    const message = `Invalid ${field} value.`
    test.each([
      ['1', true],
      ['1.0', true],
      ['-1.0', true],
      ['-1.0e5', true],
      ['-1.0e-5', true],
      ['hello', false],
      [' ', false]
    ])(
      `field: ${field}, input: %s, valid: %s`,
      async (input, valid) => {
        renderSearchEntry(<InputRange visible quantity={quantity} disableHistogram={!histogram}/>)
        const user = userEvent.setup()
        const field = await screen.findByDisplayValue(value)
        await user.clear(field)
        await user.type(field, input)
        await user.keyboard('[Enter]')
        if (valid) {
          expect(screen.queryByText(message)).toBeNull()
        } else {
          expect(screen.queryByText(message)).toBeInTheDocument()
        }
      }
    )
  }
})

describe('test histograms with only one value', () => {
  beforeAll(async () => {
    await startAPI(
      'tests.states.search.histograms_one_value',
      'tests/data/search/inputrange-one-value'
    )
  })
  afterAll(() => closeAPI())

  test.each([
    ['results.material.n_elements', 1],
    ['results.properties.electronic.band_structure_electronic.band_gap.value', 0.5],
    ['upload_create_time', 1585872000000]
  ])('quantity: %s', async (quantity, value) => {
    renderSearchEntry(<InputRange visible quantity={quantity} disableHistogram={false}/>)
    const data = defaultFilterData[quantity]
    const dtype = data.dtype

    // Check that both text fields show the only available value
    const inputValue = dtype === DType.Timestamp
      ? format(value, 'dd/MM/yyyy kk:mm')
      : value
    const inputs = await screen.findAllByDisplayValue(inputValue)
    expect(inputs.length).toBe(2)

    // Check that slider is disabled: trying to modify the sliders does not
    // update the input fields.
    const sliders = screen.getAllByRole('slider')
    const sliderMin = sliders[0]
    const sliderMax = sliders[1]
    fireEvent.keyDown(sliderMin, {key: 'Up', code: 'Up'})
    fireEvent.keyDown(sliderMax, {key: 'Down', code: 'Down'})
    const inputsNew = await screen.findAllByDisplayValue(inputValue)
    expect(inputsNew.length).toBe(2)
  })
})

test.each([
  discrete,
  discrete_histogram,
  continuous,
  continuous_histogram
])('inputs react to slider change: quantity: %s, histogram: %s', async (quantity, histogram, min, max, step) => {
  await startAPI(stateName, `tests/data/search/inputrange-${quantity}-${histogram}-slider-change`)
  renderSearchEntry(<InputRange visible quantity={quantity} disableHistogram={!histogram}/>)
  const inputMin = await screen.findByDisplayValue(min)
  const inputMax = await screen.findByDisplayValue(max)
  const sliders = screen.getAllByRole('slider')
  const sliderMin = sliders[0]
  const sliderMax = sliders[1]

  // Moving min slider changes min field
  fireEvent.keyDown(sliderMin, {key: 'Up', code: 'Up'})
  expect(inputMin.value).toBe(formatNumber(min + step))

  // Moving max slider changes max field
  fireEvent.keyDown(sliderMax, {key: 'Down', code: 'Down'})
  expect(inputMax.value).toBe(formatNumber(max - step))
  closeAPI()
})

test.each([
  continuous,
  continuous_histogram,
  discrete,
  discrete_histogram
])('sliders react to input field change: quantity: %s, histogram: %s', async (quantity, histogram, min, max) => {
  await startAPI(stateName, `tests/data/search/inputrange-${quantity}-${histogram}-input-change`)
  renderSearchEntry(<InputRange visible quantity={quantity} disableHistogram={!histogram}/>)
  const inputMin = await screen.findByDisplayValue(min)
  const inputMax = await screen.findByDisplayValue(max)
  const sliders = screen.getAllByRole('slider')
  const sliderMin = sliders[0]
  const sliderMax = sliders[1]

  // Initially both sliders at the ends
  expect(sliderMin).toHaveStyle(`left: 0%`)
  expect(sliderMax).toHaveStyle(`left: 100%`)

  // After changing the min field to 25%, min slider moves to 25% position.
  await testSliderMove(quantity, min, max, inputMin, sliderMin, 25, false, histogram)

  // After changing the max field to 75%, min slider moves to 75% position.
  await testSliderMove(quantity, min, max, inputMax, sliderMax, 75, true, histogram)
  closeAPI()
})

/**
 * Tests that a slider moves to the given location when text input changes.
 * @param {string} quantity The quantity name
 * @param {number} min Minimum value of the slider
 * @param {number} max Maximum value of the slider
 * @param {*} input Text input element
 * @param {*} slider MUI slider knob element
 * @param {number} percentage The percentage to move to.
 * @param {bool} isMax Is the max knob being moved.
 * @param {bool} isMax Is the slider shown for a histogram.
 */
async function testSliderMove(quantity, min, max, input, slider, percentage, isMax, histogram) {
  const data = defaultFilterData[quantity]
  const dtype = data.dtype
  const discretization = (histogram && dtype === DType.Int) ? 1 : 0
  const range = max - min + discretization
  const value = min + range * (percentage / 100) - (isMax ? discretization : 0)

  const user = userEvent.setup()
  await user.clear(input)
  await user.type(input, value.toString())
  await user.keyboard('[Enter]')

  const style = window.getComputedStyle(slider)
  const left = parseFloat(style.getPropertyValue('left').slice(0, -1))
  expect(left).toBeCloseTo(percentage, 8)
}
