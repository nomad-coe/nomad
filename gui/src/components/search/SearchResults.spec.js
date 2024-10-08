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

import React, {useMemo} from 'react'
import { render, screen } from '../conftest.spec'
import { SearchResults } from './SearchResults'
import { SearchContextRaw } from './SearchContext'
import { Filter } from './Filter'
import { DType } from '../../utils'

// Use a mocked SearchContext
const mockUseMemo = useMemo
jest.mock('./SearchContext', () => ({
    ...jest.requireActual('./SearchContext'),
    useSearchContext: () => ({
      ...jest.requireActual('./SearchContext').useSearchContext(),
      useResults: jest.fn(() => {
        const response = mockUseMemo(() => {
          return {
            data: [{
              my_boolean: true,
              my_float: 1.2345678,
              my_float_array: [1.2345, 2.3456],
              repeated_section: [
                {my_float: 1.2345},
                {my_float: 2.3456}
              ],
              repeated_section_outer: [
                {
                  repeated_section_inner: [
                    {my_float: 1.2345},
                    {my_float: 2.3456}
                  ]
                },
                {
                  repeated_section_inner: [
                    {my_float: 3.456}
                  ]
                }
              ],
              data: {
                my_float: 1.2345678
              }
            }],
            pagination: {total: 1},
            setPagination: () => {}
          }
        }, [])
        return response
      })
    })
}))

describe('test search columns', () => {
  test.each([
    ['boolean scalar', {quantity: 'my_boolean', selected: true}, 'My boolean', 'True', true],
    ['float scalar', {quantity: 'my_float', selected: true, format: {decimals: 2}}, 'My float', '1.23', true],
    ['float from repeating section', {quantity: 'repeated_section[*].my_float', selected: true, format: {decimals: 2}}, 'My float', '1.23, 2.35', false],
    ['float from nested repeating section', {quantity: 'repeated_section_outer[*].repeated_section_inner[*].my_float', selected: true, format: {decimals: 2}}, 'My float', '1.23, 2.35, 3.46', false],
    ['float from array', {quantity: 'my_float_array', selected: true, format: {decimals: 2}}, 'My float array', '1.23, 2.35', false],
    ['custom title', {quantity: 'my_float', title: 'Custom', selected: true, format: {decimals: 2}}, 'Custom', '1.23', true],
    ['custom quantity', {quantity: 'data.my_float#my_plugin.schema_packages.my_schema.MySchema.my_float', selected: true, format: {decimals: 2}}, 'My float', '1.23', true]
  ])('%s', async (name, column, expected_title, expected_value, sortable) => {
    render(
      <SearchContextRaw
        resource="entries"
        id='entries'
        initialFilterData={{
          'my_float': new Filter(undefined, {quantity: 'my_float', dtype: DType.Float}),
          'my_boolean': new Filter(undefined, {quantity: 'my_boolean', dtype: DType.Boolean}),
          'my_float_array': new Filter(undefined, {quantity: 'my_float_array', dtype: DType.Float, shape: '*'}),
          'repeated_section.my_float': new Filter(undefined, {quantity: 'repeated_section.my_float', dtype: DType.Float}),
          'repeated_section_outer.repeated_section_inner.my_float': new Filter(undefined, {quantity: 'repeated_section_outer.repeated_section_inner.my_float', dtype: DType.Float}),
          'data.my_float#my_plugin.schema_packages.my_schema.MySchema.my_float': new Filter(undefined, {quantity: 'data.my_float#my_plugin.schema_packages.my_schema.MySchema.my_float', dtype: DType.Float})
        }}
        initialColumns={[column]}
      >
        <SearchResults />
      </SearchContextRaw>
    )
    screen.getByText(expected_title)
    screen.getByText(expected_value)
    const sortButton = screen.queryByTestId(`sortable_${column.quantity}`)
    if (sortable) {
      expect(sortButton).not.toBeNull()
    } else {
      expect(sortButton).toBeNull()
    }
  })
})
