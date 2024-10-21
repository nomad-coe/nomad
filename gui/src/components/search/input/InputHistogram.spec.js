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
import React, { useCallback, useMemo, useState } from 'react'
import { fireEvent } from '@testing-library/dom'
import { format } from 'date-fns'
import { render, screen } from '../../conftest.spec'
import { Filter } from '../Filter'
import InputHistogram from './InputHistogram'
import userEvent from '@testing-library/user-event'
import { DType, formatNumber } from '../../../utils'
import { useSearchContext, SearchContextRaw } from '../SearchContext'

// We set an initial mock for the SearchContext module
const mockSetFilter = jest.fn()
jest.mock('../SearchContext', () => ({
    ...jest.requireActual('../SearchContext'),
    useSearchContext: jest.fn()
}))

describe('', () => {
  // Provide a default implementation for the search context mock
  beforeEach(() => {
    useSearchContext.mockImplementation(() => ({
        ...jest.requireActual('../SearchContext').useSearchContext(),
        useAgg: (quantity, visible, id, config) => {
          return useMemo(() => {
            return visible
              ? {
                  data: [{value: 0, count: 10}, {value: 1, count: 20}, {value: 2, count: 20}],
                  interval: 1
                }
              : undefined
          }, [visible])
        },
        useFilterState: (quantity) => {
          const response = useMemo(() => {
            return [undefined, mockSetFilter]
          }, [])
          return response
        }
    }))
  })

  describe('test showHeader', () => {
    test.each([
      ['show header', {showHeader: true}],
      ['do not show header', {showHeader: false}]
    ])('%s', async (name, config) => {
      renderInputHistogram(config, new Filter(undefined, {quantity: 'test'}))
      if (config.showHeader) {
        expect(screen.getByText('Test')).toBeInTheDocument()
      } else {
        expect(screen.queryByText('Test')).not.toBeInTheDocument()
      }
    })
  })

  describe('test title', () => {
    test.each([
      ['default title', {}, 'Test'],
      ['custom title', {title: 'Custom title'}, 'Custom title']
    ])('%s', async (name, config, expected) => {
      renderInputHistogram(config, new Filter(undefined, {quantity: 'test'}))
      expect(screen.getByText(expected)).toBeInTheDocument()
    })
  })

  describe('test showStatistics', () => {
    test.each([
      ['show statistics', {showStatistics: true}],
      ['do not show statistics', {showStatistics: false}]
    ])('%s', async (name, config) => {
      renderInputHistogram(config, new Filter(undefined, {quantity: 'test'}))
    const option = screen.queryAllByText('1')
    const scaling = screen.queryByText('linear')
    if (config.showStatistics) {
      expect(option).toHaveLength(2)
      expect(scaling).toBeInTheDocument()
    } else {
      expect(option).toHaveLength(0)
      expect(scaling).not.toBeInTheDocument()
    }
    })
  })

  describe('test showInput', () => {
    test.each([
      ['show input', {showInput: true}],
      ['do not show input', {showInput: false}]
    ])('%s', async (name, config) => {
      renderInputHistogram(config, new Filter(undefined, {quantity: 'test'}))
      if (config.showInput) {
        expect(screen.getByText('min:')).toBeInTheDocument()
        expect(screen.getByText('max:')).toBeInTheDocument()
      } else {
        expect(screen.queryByPlaceholderText('min:')).not.toBeInTheDocument()
        expect(screen.queryByPlaceholderText('max:')).not.toBeInTheDocument()
      }
    })
  })

  describe('test invalid/valid numeric input', () => {
    for (const isMin of [false, true]) {
      const value = isMin ? 0 : 3
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
        `input: %s, valid: %s`,
        async (input, valid) => {
          const user = userEvent.setup()
          renderInputHistogram({}, new Filter(undefined, {quantity: 'test'}))
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

  describe('test input field change', () => {
    // Custom aggregation data for this test
    beforeEach(() => {
      useSearchContext.mockImplementation(() => ({
          ...jest.requireActual('../SearchContext').useSearchContext(),
          useAgg: (quantity, visible, id, config) => {
            return useMemo(() => {
              return visible
                ? {
                    data: [{value: 0, count: 10}, {value: 1, count: 20}, {value: 2, count: 20}],
                    interval: 1
                  }
                : undefined
            }, [visible])
          },
          useFilterState: (quantity) => {
            const [filter, setFilter] = useState()
            const setFilterWrapper = useCallback((value) => {
              setFilter(value)
              mockSetFilter(value)
            }, [])
            const response = useMemo(() => {
              return [filter, setFilterWrapper]
            }, [filter, setFilterWrapper])
            return response
          }
      }))
    })
    test.each([
      ['min field', 0, 1, 3, 0, `left: 0%`],
      ['max field', 3, 0, 1, 1, `left: 100%`]
    ])('%s', async (name, value, gte, lte, slider, position) => {
      const filter = new Filter(undefined, {quantity: 'test'})
      const user = userEvent.setup()
      renderInputHistogram({}, filter)
      const sliders = screen.getAllByRole('slider')
      slider = sliders[slider]

      // Initially slider at the end
      expect(slider).toHaveStyle(position)

      // Type in a new value
      const input = await screen.findByDisplayValue(value)
      await user.clear(input)
      await user.type(input, '1')
      await user.keyboard('[Enter]')

      // setFilter is triggered
      expect(mockSetFilter.mock.calls).toHaveLength(1)
      const argument = mockSetFilter.mock.calls[0][0]()
      expect(argument.gte.value()).toBe(gte)
      expect(argument.lte.value()).toBe(lte)

      // Changing input moves slider
      await testSliderMove(filter.dtype, 0, 3, input, slider, 50, false)
    })
  })

  describe('test slider change', () => {
    test.each([
      ['min slider', 0, 1, 3, 0, 1, {key: 'Up', code: 'Up'}],
      ['max slider', 3, 0, 2, 1, -1, {key: 'Down', code: 'Down'}]
    ])('%s', async (name, value, gte, lte, slider, step, input) => {
      renderInputHistogram({}, new Filter(undefined, {quantity: 'test'}))
      const sliders = screen.getAllByRole('slider')
      slider = sliders[slider]

      // Change slider value with arrow keys
      const inputMin = await screen.findByDisplayValue(value)
      fireEvent.keyDown(slider, input)

      // setFilter is triggered
      expect(mockSetFilter.mock.calls).toHaveLength(1)
      const argument = mockSetFilter.mock.calls[0][0]
      expect(argument.gte.value()).toBe(gte)
      expect(argument.lte.value()).toBe(lte)

      // Moving slider changes input field
      expect(mockSetFilter.mock.calls).toHaveLength(1)
      expect(inputMin.value).toBe(formatNumber(value + step))
    })
  })

  describe('test histograms with only one value', () => {
    // Custom aggregation data for this test
    beforeEach(() => {
      useSearchContext.mockImplementation(() => ({
          ...jest.requireActual('../SearchContext').useSearchContext(),
          useAgg: (quantity, visible, id, config) => {
              return useMemo(() => {
                  return visible
                      ? { data: [{ value: 1, count: 10 }], interval: 0 }
                      : undefined
              }, [visible])
          }
      }))
    })
    test.each([
      ['integer', new Filter(undefined, {quantity: 'test', dtype: DType.Int}), 1],
      ['float', new Filter(undefined, {quantity: 'test', dtype: DType.Float}), 1],
      ['timestamp', new Filter(undefined, {quantity: 'test', dtype: DType.Timestamp}), 1]
    ])('quantity: %s', async (name, filter, value) => {
      renderInputHistogram({}, filter)

      // Check that both text fields show the only available value
      const inputValue = filter.dtype === DType.Timestamp
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
async function testSliderMove(dtype, min, max, input, slider, percentage, isMax) {
  const discretization = (dtype === DType.Int) ? 1 : 0
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

// Helper function for rendering
function renderInputHistogram(config, filter) {
  const searchQuantities = {test: filter}
  render(
    <SearchContextRaw
      resource="entries"
      id='entries'
      initialSearchQuantities={searchQuantities}
    >
      <InputHistogram visible x={{search_quantity: filter.quantity}} {...config}/>
    </SearchContextRaw>
  )
}
