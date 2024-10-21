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
import React, { useMemo } from 'react'
import { waitFor, within } from '@testing-library/dom'
import { render, screen } from '../../conftest.spec'
import { expectInputItem } from '../conftest.spec'
import { SearchContextRaw } from '../SearchContext'
import { Filter } from '../Filter'
import { isNumber } from 'lodash'
import InputTerms from './InputTerms'
import userEvent from '@testing-library/user-event'

// Use a mocked SearchContext
const mockSetFilter = jest.fn()
const mockUseMemo = useMemo
jest.mock('../SearchContext', () => ({
    ...jest.requireActual('../SearchContext'),
    useSearchContext: () => ({
      ...jest.requireActual('../SearchContext').useSearchContext(),
      useAgg: (quantity, visible, id, config) => {
        const response = mockUseMemo(() => {
          return visible
            ? {data: [
              {value: 'A', count: 6},
              {value: 'B', count: 5},
              {value: 'C', count: 4},
              {value: 'D', count: 3},
              {value: 'E', count: 2},
              {value: 'F', count: 1}
            ].slice(0, config.size)}
            : undefined
        }, [])
        return response
      },
      useFilterState: jest.fn((quantity) => {
        const response = mockUseMemo(() => {
          return [undefined, mockSetFilter]
        }, [])
        return response
      })
    })
}))

describe('test options', () => {
  test.each([
    [
      'show all options for enum by default',
      {options: undefined},
      new Filter(undefined, {quantity: 'test', options: {A: {label: 'A'}, B: {label: 'B'}}}),
      [{label: 'A'}, {label: 'B'}]
    ],
    [
      'show 5 options for str by default',
      {options: undefined},
      new Filter(undefined, {quantity: 'test'}),
      [{label: 'A'}, {label: 'B'}, {label: 'C'}, {label: 'D'}, {label: 'E'}]
    ],
    [
      'no options',
      {options: 0},
      new Filter(undefined, {quantity: 'test'}),
      []
    ],
    [
      'limited options',
      {options: 2},
      new Filter(undefined, {quantity: 'test'}),
      [{label: 'A'}, {label: 'B'}]
    ],
    [
      'custom options',
      {options: {B: {label: 'B'}}},
      new Filter(undefined, {quantity: 'test', options: {A: {label: 'A'}, B: {label: 'B'}}}),
      [{label: 'B'}]
    ]
  ])('%s', async (name, config, filter, expected) => {
    renderInputTerms(config, filter)

    // Wait for possible placeholder to disappear
    await waitFor(() => expect(screen.queryByTestId(`input-terms-placeholder`)).toBe(null))

    // Check that each expected item appears
    expect(screen.queryAllByRole('checkbox').length).toBe(expected.length)
    for (const item of expected) {
      expectInputItem(item)
    }

    // When getting options dynamically, test that the "show more" button is
    // shown, but "show less" is not shown
    if (isNumber(config.options) && config.options > 0) {
      expect(screen.getByText('Show more')).toBeInTheDocument()
      expect(screen.queryByText('Show less')).not.toBeInTheDocument()
    }
  })
})

describe('test showHeader', () => {
  test.each([
    ['show header', {showHeader: true}],
    ['do not show header', {showHeader: false}]
  ])('%s', async (name, config) => {
    renderInputTerms(config, new Filter(undefined, {quantity: 'test'}))
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
    renderInputTerms(config, new Filter(undefined, {quantity: 'test'}))
    expect(screen.getByText(expected)).toBeInTheDocument()
  })
})

describe('test showStatistics', () => {
  test.each([
    ['show statistics', {showStatistics: true}],
    ['do not show statistics', {showStatistics: false}]
  ])('%s', async (name, config) => {
    renderInputTerms(config, new Filter(undefined, {quantity: 'test'}))
    const option = screen.queryByText('6')
    const scaling = screen.queryByText('linear')
    if (config.showStatistics) {
      expect(option).toBeInTheDocument()
      expect(scaling).toBeInTheDocument()
    } else {
      expect(option).not.toBeInTheDocument()
      expect(scaling).not.toBeInTheDocument()
    }
  })
})

describe('test showInput', () => {
  test.each([
    ['show input', {showInput: true}],
    ['do not show input', {showInput: false}]
  ])('%s', async (name, config) => {
    renderInputTerms(config, new Filter(undefined, {quantity: 'test'}))
    if (config.showInput) {
      expect(screen.getByPlaceholderText('Type here')).toBeInTheDocument()
    } else {
      expect(screen.queryByPlaceholderText('Type here')).not.toBeInTheDocument()
    }
  })
})

test.only('test item selection', async () => {
  renderInputTerms({}, new Filter(undefined, {quantity: 'test'}))

  // Wait for possible placeholder to disappear
  await waitFor(() => expect(screen.queryByTestId(`input-terms-placeholder`)).toBe(null))

  // Select one item
  const checkbox = queryByInputItemName('A')
  await userEvent.click(checkbox)

  // Check that the setFilter function is called once with the correct argument
  expect(mockSetFilter.mock.calls).toHaveLength(1)
  const argument = mockSetFilter.mock.calls[0][0]
  const expectedArgument = new Set(['A'])
  expect(argument.size === expectedArgument.size).toBe(true)
  expect([...argument].every((x) => expectedArgument.has(x))).toBe(true)
})

// Helper function for rendering
function renderInputTerms(config, filter) {
  const searchQuantities = {test: filter}
  render(
    <SearchContextRaw
      resource="entries"
      id='entries'
      initialSearchQuantities={searchQuantities}
    >
      <InputTerms visible searchQuantity={filter.quantity} {...config}/>
    </SearchContextRaw>
  )
}

/**
 * Finds the checkbox corresponding to an InputItem with the given value.
 * @param {string} name The option value that is displayed
 * @returns {element} The checkbox input HTML element.
 */
function queryByInputItemName(option, root = screen) {
  const inputLabel = root.queryByText(option)
  const inputCheckbox = inputLabel && within(inputLabel.closest('label')).getByRole('checkbox')
  return inputCheckbox
}
