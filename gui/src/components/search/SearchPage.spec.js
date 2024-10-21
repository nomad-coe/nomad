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
import { expectMenu, expectSearchResults } from './conftest.spec'
import { ui } from '../../config'
import { SearchContext } from './SearchContext'
import SearchPage from './SearchPage'

// We set an initial mock for the SearchContext module
const mockSetFilter = jest.fn()
const mockUseMemo = useMemo
jest.mock('./SearchContext', () => ({
    ...jest.requireActual('./SearchContext'),
    useSearchContext: () => ({
      ...jest.requireActual('./SearchContext').useSearchContext(),
      useAgg: (quantity, visible, id, config) => {
        const response = mockUseMemo(() => {
          return undefined
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
          return {data: [{}], pagination: {total: 1}}
        }, [])
        return response
      })
    })
}))

describe('', () => {
  test('render search page components', async () => {
    const context = ui.apps.options.entries
    render(
      <SearchContext
          resource={context.resource}
          initialPagination={context.pagination}
          initialColumns={context.columns}
          initialRows={context.rows}
          initialMenu={context.menu}
          initialFiltersLocked={context.filters_locked}
          initialDashboard={context?.dashboard}
          initialSearchSyntaxes={context?.search_syntaxes}
          id={context?.path}
      >
        <SearchPage />
      </SearchContext>
    )
    // Test that menu is shown
    await expectMenu(context.menu)

    // Test that search bar is shown
    screen.getByPlaceholderText('Type your query or keyword here')

    // Test that query is shown
    screen.getByText('Your query will be shown here')

    // Test that results table is shown
    await expectSearchResults(context.columns)
  })
})
