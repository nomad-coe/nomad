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
import { renderHook, act } from '@testing-library/react-hooks'
import { WrapperSearch } from './conftest.spec'
import { useSearchContext } from './SearchContext'
import { Quantity } from '../units/Quantity'
import { isEqualWith } from 'lodash'
import { SearchSuggestion, SuggestionType } from './SearchSuggestion'
import { WrapperDefault } from '../conftest.spec'

describe('parseQuery', function() {
  test.each([
    ['unit not specified', 'results.material.topology.cell.a', '1', new Quantity(1, 'angstrom'), undefined],
    ['unit specified', 'results.material.topology.cell.a', '1 m', new Quantity(1, 'meter'), undefined],
    ['cannot parse number', 'results.material.topology.cell.a', 'a', undefined, 'Could not parse the number a'],
    ['filter hat accepts multiple values is wrapped in set', 'results.material.material_id', 'abcd', new Set(['abcd']), undefined],
    ['filter that does not accept multiple values is not wrapped in set', 'visibility', 'public', 'public', undefined]
  ])('%s', async (name, quantity, input, output, error) => {
    const { result: resultUseSearchContext } = renderHook(() => useSearchContext(), { wrapper: WrapperSearch })
    const { result: resultUseParseQuery } = renderHook(() => resultUseSearchContext.current.useParseQuery(), {})
    const parseQuery = resultUseParseQuery.current
    if (!error) {
      function customizer(a, b) {
        if (a instanceof Quantity) {
          return a.equal(b)
        }
      }
      expect(isEqualWith(parseQuery(quantity, input), output, customizer)).toBe(true)
    } else {
      expect(() => parseQuery(quantity, input)).toThrow(error)
    }
    }
  )
})

const suggestionsInitial = [
  new SearchSuggestion({input: 'old', type: SuggestionType.Freetext})
]
describe('suggestions', function() {
  beforeAll(() => {
    window.localStorage.setItem(
      'nomad-searchcontext-entries',
      JSON.stringify(suggestionsInitial)
    )
  })
  afterAll(() => {
    window.localStorage.removeItem('nomad-searchcontext-entries')
  })
  test('load initial suggestions from localStorage', async () => {
    const { result: resultUseSearchContext } = renderHook(() => useSearchContext(), { wrapper: WrapperSearch })
    const { result: resultUseSearchSuggestions } = renderHook(() => resultUseSearchContext.current.useSearchSuggestions(''), { wrapper: WrapperSearch })
    const suggestions = resultUseSearchSuggestions.current
    expect(suggestions).toHaveLength(suggestionsInitial.length)
    suggestions.forEach((value, index) => {
      expect(value.key === suggestionsInitial[index].key)
    })
  })
  test('save new suggestion to localStorage', async () => {
    const { result: resultUseSearchContext } = renderHook(() => useSearchContext(), { wrapper: WrapperSearch })
    const { result: resultusePushSearchSuggestion } = renderHook(() => resultUseSearchContext.current.usePushSearchSuggestion(), { wrapper: WrapperSearch })
    const pushSearchSuggestion = resultusePushSearchSuggestion.current
    const newSuggestion = new SearchSuggestion({input: 'new', type: SuggestionType.Freetext})
    act(() => {
      pushSearchSuggestion(newSuggestion)
    })
    const suggestions = JSON.parse(window.localStorage.getItem('nomad-searchcontext-entries')).map(x => new SearchSuggestion(x))
    expect(suggestions).toHaveLength(suggestionsInitial.length + 1)
    expect(suggestions[0].key).toBe(newSuggestion.key)
  })
})

describe('reading query from URL', function() {
  test.each([
    ['no query parameters', '', {}],
    ['query parameter', '?upload_id=testing', {upload_id: new Set(['testing'])}],
    ['query parameter with keycloak iss', '?upload_id=testing#iss=https://nomad-lab.eu/fairdi/keycloak/auth/realms/fairdi_nomad_prod', {upload_id: new Set(['testing'])}],
    ['query parameter with keycloak error', '?upload_id=testing#error=login_required', {upload_id: new Set(['testing'])}],
    ['query parameter with keycloak state', '?upload_id=testing#state=here is a state', {upload_id: new Set(['testing'])}]
  ])('%s', async (name, params, expected_query) => {
    // Set window.location.href
    // eslint-disable-next-line no-global-assign
    window = Object.create(window)
    const url = "http://testing.com" + params
    Object.defineProperty(window, 'location', {
      value: {href: url},
      writable: true
    })

    // Call hooks to check that query is correctly read
    const { result: resultUseSearchContext } = renderHook(() => useSearchContext(), { wrapper: WrapperSearch })
    const { result: resultUseQuery } = renderHook(() => resultUseSearchContext.current.useQuery(), { wrapper: WrapperDefault})
    const query = resultUseQuery.current
    expect(query).toMatchObject(expected_query)
  })
})
