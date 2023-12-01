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
import { renderHook } from '@testing-library/react-hooks'
import { WrapperSearch } from './conftest.spec'
import { useSearchContext } from './SearchContext'
import { Quantity } from '../units/Quantity'
import { isEqualWith } from 'lodash'

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
