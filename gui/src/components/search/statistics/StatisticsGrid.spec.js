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
import { render, screen, startAPI, closeAPI } from '../../conftest.spec'
import userEvent from '@testing-library/user-event'
import StatisticsGrid from './StatisticsGrid'
import StatisticsToggle from './StatisticsToggle'
import { filterData } from '../FilterRegistry'
import { SearchContext } from '../SearchContext'
import {
  expectInputList,
  expectInputPeriodicTable,
  expectInputRange
} from '../conftest.spec'

// Mock the useResizeObserver hook. The test environment does not provide any
// resize events on it's own.
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {width: 1000, height: 300, ref: undefined}
  }}
})

describe('displaying, removing and re-adding an item in the grid works', () => {
  beforeAll(async () => {
    await startAPI('tests.states.search.search', 'tests/data/search/statisticsgrid')
  })
  afterAll(() => closeAPI())

  // Each unique component that can be inserted in the grid should be tested
  // here.
  test.each([
    [
      'InputList',
      'results.material.structural_type',
      async (quantity, loaded) => await expectInputList(quantity, loaded, ['2D', 'bulk', 'molecule / cluster'], 'all')
    ],
    [
      'InputRange',
      'results.material.n_elements',
      async (quantity, loaded) => await expectInputRange(quantity, loaded, true, true)
    ],
    [
      'InputPeriodicTable',
      'results.material.elements',
      async (quantity, loaded) => await expectInputPeriodicTable(quantity, loaded, ['H', 'C', 'N', 'Ti', 'Zr', 'Nb', 'I', 'Hf', 'Ta', 'Pb'])
    ]
  ])('%s', async (component, quantity, test) => {
    const initialStats = {[quantity]: {index: 0}}
    render(
      <SearchContext resource="entries" initialStatistics={initialStats}>
        <StatisticsToggle quantity={quantity} data-testid="statisticstoggle"/>
        <StatisticsGrid />
      </SearchContext>
    )

    // Assert that the item is displayed initially
    await test(quantity, false)

    // Remove item, check that it is gone. A test id is used to fetch the
    // statistics toggle since it is an icon where the only identifier is the
    // tooltip text which may appear in several locations depending on whether
    // the button is focused or not.
    const toggleButton = screen.getByTestId('statisticstoggle')
    const data = filterData[quantity]
    const label = data.label
    expect(screen.queryByText(label, {exact: false})).toBeInTheDocument()
    await userEvent.click(toggleButton)
    expect(screen.queryByText(label, {exact: false})).not.toBeInTheDocument()

    // Re-add item, check that it appears. This time the data is already loaded.
    await userEvent.click(toggleButton)
    await test(quantity, true)
  })
})
