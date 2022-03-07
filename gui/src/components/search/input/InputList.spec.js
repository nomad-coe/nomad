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
import { waitForElementToBeRemoved } from '@testing-library/dom'
import { startAPI, closeAPI, screen } from '../../conftest'
import { renderSearchEntry, expectInputHeader } from '../conftest'
import InputList from './InputList'

// Mock the useResizeObserver hook. The test environment does not provide any
// resize events on it's own.
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {height: 300, ref: undefined}
  }}
})

const quantity = 'results.material.structural_type'
const stateName = 'tests.states.search.search'

describe('', () => {
  beforeEach(() => {
    // API state with single terms aggregation result
    startAPI(stateName, 'tests/data/search/terms_aggregation_structural_type')

    // Render InputList within an entry search context. The component is wrapped
    // inside a div that controls the final size.
    renderSearchEntry(<InputList
      quantity={quantity}
      visible
      draggable
      aggId="statistics"
      data-testid="inputlist"
    />
    )
  })
  afterEach(() => closeAPI())

  test('initial state is loaded correctly', async () => {
    // Test immediately displayed elements
    expectInputHeader(quantity)

    // Test that placeholder is shown while loading
    const placeholder = screen.getByTestId('inputlist-placeholder')
    expect(placeholder).toBeInTheDocument()

    // Check that placeholder disappears
    await waitForElementToBeRemoved(() => screen.getByTestId('inputlist-placeholder'))

    // Test elements that are displayed after API response
    expect(await screen.findByText('molecule / cluster')).toBeInTheDocument()
    expect(await screen.findByText('2D')).toBeInTheDocument()
    expect(await screen.findByText('bulk')).toBeInTheDocument()
    expect(screen.getByText('Showing all 3 items')).toBeInTheDocument()
  })
})
