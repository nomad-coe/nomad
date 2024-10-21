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
import { render, screen } from '../../conftest.spec'
import { Filter } from '../Filter'
import InputPeriodicTable from './InputPeriodicTable'
import userEvent from '@testing-library/user-event'
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
              ? {data: [{value: 'H', count: 123}]}
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
      renderPeriodicTable(config, new Filter(undefined, {quantity: 'test'}))
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
      renderPeriodicTable(config, new Filter(undefined, {quantity: 'test'}))
      expect(screen.getByText(expected)).toBeInTheDocument()
    })
  })

  describe('test showStatistics', () => {
    test.each([
      ['show statistics', {showStatistics: true}],
      ['do not show statistics', {showStatistics: false}]
    ])('%s', async (name, config) => {
      renderPeriodicTable(config, new Filter(undefined, {quantity: 'test'}))

      const option = screen.queryAllByText('123')
      const scaling = screen.queryByText('linear')
      if (config.showStatistics) {
        expect(option).toHaveLength(1)
        expect(scaling).toBeInTheDocument()
      } else {
        expect(option).toHaveLength(0)
        expect(scaling).not.toBeInTheDocument()
      }
    })
  })

  describe('test selection', () => {
    test.each([
      ['show statistics', {showStatistics: true}],
      ['do not show statistics', {showStatistics: false}]
    ])('%s', async (name, config) => {
      renderPeriodicTable(config, new Filter(undefined, {quantity: 'test'}))

      // Click on hydrogen
      const cButton = screen.getByTestId('Hydrogen')
      await userEvent.click(cButton)

      // setFilter is called
      expect(mockSetFilter.mock.calls).toHaveLength(1)
      const argument = mockSetFilter.mock.calls[0][0]()
      const expectedArgument = new Set(['H'])
      expect(argument.size === expectedArgument.size).toBe(true)
      expect([...argument].every((x) => expectedArgument.has(x))).toBe(true)
    })
  })
})

// Helper function for rendering
function renderPeriodicTable(config, filter) {
  const searchQuantities = {test: filter}
  render(
    <SearchContextRaw
      resource="entries"
      id='entries'
      initialSearchQuantities={searchQuantities}
    >
      <InputPeriodicTable visible searchQuantity={filter.quantity} {...config}/>
    </SearchContextRaw>
  )
}
