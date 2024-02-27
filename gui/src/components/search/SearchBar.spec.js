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
import { ui } from '../../config'
import { render, screen, within } from '../conftest.spec'
import userEvent from '@testing-library/user-event'
import { SearchContext } from './SearchContext'
import SearchBar from './SearchBar'
import { SearchSuggestion, SuggestionType } from './SearchSuggestion'

describe('searchbar queries', function() {
  test.each([
    ['equality', 'results.material.n_elements = 1', 'equality'],
    ['existence', 'results.material.n_elements = *', 'existence'],
    ['range, -gt', 'results.material.n_elements > 0', 'range'],
    ['range, -gte', 'results.material.n_elements >= 0', 'range'],
    ['range, -lt', 'results.material.n_elements < 0', 'range'],
    ['range, -lte', 'results.material.n_elements <= 0', 'range'],
    ['range, gt-', '0 > results.material.n_elements', 'range'],
    ['range, gte-', '0 >= results.material.n_elements', 'range'],
    ['range, lt-', '0 < results.material.n_elements', 'range'],
    ['range, lte', '0 <= results.material.n_elements', 'range'],
    ['range, lt-lt', '0 < results.material.n_elements < 1', 'range'],
    ['range, lte-lt', '0 <= results.material.n_elements < 1', 'range'],
    ['range, lt-lte', '0 < results.material.n_elements <= 1', 'range'],
    ['range, lte-lte', '0 <= results.material.n_elements <= 1', 'range'],
    ['range, gt-gt', '0 > results.material.n_elements > 1', 'range'],
    ['range, gte-gt', '0 >= results.material.n_elements > 1', 'range'],
    ['range, gt-gte', '0 > results.material.n_elements >= 1', 'range'],
    ['range, gte-gte', '0 >= results.material.n_elements >= 1', 'range']
  ])('%s', async (name, input, type) => {
    const context = ui.apps.options.entries
    render(
      <SearchContext
          resource={context.resource}
          initialPagination={context.pagination}
          initialColumns={context.columns}
          initialRows={context.rows}
          initialFilterMenus={context.filter_menus}
          initialFiltersLocked={context.filters_locked}
          initialDashboard={context?.dashboard}
          initialSearchSyntaxes={context?.search_syntaxes}
      >
        <SearchBar />
      </SearchContext>
    )
    const textInput = screen.getByRole('textbox')
    await userEvent.type(textInput, input)
    expect(screen.getByRole('textbox')).toHaveValue(input)
    await userEvent.keyboard('{enter}')
    expect(screen.getByRole('textbox')).toHaveValue('')
  })
})

const suggestionsInitial = [
  new SearchSuggestion({input: 'old', type: SuggestionType.Freetext, history: true}),
  new SearchSuggestion({input: 'older', type: SuggestionType.Freetext, history: true}),
  new SearchSuggestion({input: 'oldest', type: SuggestionType.Freetext, history: true})
]
describe('suggestions: history', function() {
  beforeAll(() => {
    window.localStorage.removeItem('nomad-searchcontext-entries')
    window.localStorage.setItem(
      'nomad-searchcontext-entries',
      JSON.stringify(suggestionsInitial)
    )
  })
  afterAll(() => {
    window.localStorage.removeItem('nomad-searchcontext-entries')
  })
  test('initial suggestions are shown in correct order upon focus', async () => {
    const context = ui.apps.options.entries
    render(
      <SearchContext
          resource={context.resource}
          initialPagination={context.pagination}
          initialColumns={context.columns}
          initialRows={context.rows}
          initialFilterMenus={context.filter_menus}
          initialFiltersLocked={context.filters_locked}
          initialDashboard={context?.dashboard}
          initialSearchSyntaxes={context?.search_syntaxes}
          id={context?.path}
      >
        <SearchBar />
      </SearchContext>
    )
    const textInput = screen.getByRole('textbox')
    await userEvent.click(textInput)
    const options = screen.getAllByRole('option')
    expect(options.length).toBe(3)
    expect(options[0]).toHaveTextContent('old')
    expect(options[1]).toHaveTextContent('older')
    expect(options[2]).toHaveTextContent('oldest')
  })
  test('options are filtered according to input', async () => {
    const context = ui.apps.options.entries
    render(
      <SearchContext
          resource={context.resource}
          initialPagination={context.pagination}
          initialColumns={context.columns}
          initialRows={context.rows}
          initialFilterMenus={context.filter_menus}
          initialFiltersLocked={context.filters_locked}
          initialDashboard={context?.dashboard}
          initialSearchSyntaxes={context?.search_syntaxes}
          id={context?.path}
      >
        <SearchBar />
      </SearchContext>
    )
    const textInput = screen.getByRole('textbox')
    await userEvent.type(textInput, 'olde')
    const options = screen.getAllByRole('option')
    expect(options.length).toBe(2)
    expect(options[0]).toHaveTextContent('older')
    expect(options[1]).toHaveTextContent('oldest')
  })
  test('clicking delete icon removes option from list', async () => {
    const context = ui.apps.options.entries
    render(
      <SearchContext
          resource={context.resource}
          initialPagination={context.pagination}
          initialColumns={context.columns}
          initialRows={context.rows}
          initialFilterMenus={context.filter_menus}
          initialFiltersLocked={context.filters_locked}
          initialDashboard={context?.dashboard}
          initialSearchSyntaxes={context?.search_syntaxes}
          id={context?.path}
      >
        <SearchBar />
      </SearchContext>
    )
    const textInput = screen.getByRole('textbox')
    await userEvent.click(textInput)
    const options = screen.getAllByRole('option')
    await userEvent.hover(options[0])
    const removeButton = within(options[0]).getByTooltip('Remove from search history')
    await userEvent.click(removeButton)
    const optionsNew = screen.getAllByRole('option')
    expect(optionsNew.length).toBe(2)
    expect(optionsNew[0]).toHaveTextContent('older')
    expect(optionsNew[1]).toHaveTextContent('oldest')
  })
})
