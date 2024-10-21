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
import { render, screen } from '../conftest.spec'
import { SearchContextRaw } from './SearchContext'
import SearchMenu from './SearchMenu'
import { Filter } from './Filter'
import { DType } from '../../utils'

// We set an initial mock for the SearchContext module
const mockSetFilter = jest.fn()
const mockUseMemo = useMemo
jest.mock('./SearchContext', () => ({
    ...jest.requireActual('./SearchContext'),
    useSearchContext: () => ({
      ...jest.requireActual('./SearchContext').useSearchContext(),
      useAgg: (quantity, visible, id, config) => {
        const response = mockUseMemo(() => {
          return visible
            ? {data: quantity === 'my_histogram'
              ? [{value: 1, count: 5}]
              : [
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
      }),
      useResults: jest.fn((quantity) => {
        const response = mockUseMemo(() => {
          return undefined
        }, [])
        return response
      })
    })
}))

describe('test menu items', () => {
  test.each([
    ['terms', {type: 'terms', search_quantity: 'my_terms'}, 'My terms'],
    ['histogram', {type: 'histogram', x: {search_quantity: 'my_histogram'}}, 'My histogram'],
    ['periodic table', {type: 'periodic_table', search_quantity: 'my_periodic_table'}, 'My periodic table'],
    ['nested object', {type: 'nested_object', path: 'my_nested_object'}, 'My nested object'],
    ['visibility', {type: 'visibility'}, 'Visibility'],
    ['definitions', {type: 'definitions'}, 'Definitions'],
    ['optimade', {type: 'optimade'}, 'Optimade query'],
    ['custom quantities', {type: 'custom_quantities'}, 'Custom quantities'],
    ['menu', {type: 'menu', 'title': 'Submenu'}, 'Submenu']
  ])('%s', async (name, menuItem, text) => {
    const mainmenu = {
      items: [menuItem]
    }
    render(
      <SearchContextRaw
        resource="entries"
        id='entries'
        initialSearchQuantities={{
          my_terms: new Filter(undefined, {quantity: 'my_terms', dtype: DType.String}),
          my_histogram: new Filter(undefined, {quantity: 'my_histogram', dtype: DType.Int}),
          my_periodic_table: new Filter(undefined, {quantity: 'my_periodic_table', dtype: DType.String}),
          my_nested_object: new Filter(undefined, {quantity: 'my_nested_object', dtype: DType.String}),
          visibility: new Filter(undefined, {quantity: 'visibility', dtype: DType.String, global: true}),
          quantities: new Filter(undefined, {quantity: 'quantities', dtype: DType.String})
        }}
        initialMenu={mainmenu}
      >
        <SearchMenu/>
      </SearchContextRaw>
    )
    // Test that the menu item is displayed
    screen.getByText(text)
  })
})
