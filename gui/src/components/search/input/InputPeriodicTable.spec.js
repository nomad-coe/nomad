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
import userEvent from '@testing-library/user-event'
import { waitFor } from '@testing-library/dom'
import { startAPI, closeAPI, screen } from '../../conftest.spec'
import { renderSearchEntry, expectInputHeader } from '../conftest.spec'
import elementData from '../../../elementData.json'
import InputPeriodicTable from './InputPeriodicTable'

const quantity = 'results.material.elements'
const stateName = 'tests.states.search.search'

/**
 * Tests that only the given elements are selectable in the periodic table.
 * @param {array} elements List of chemical symbols.
 * @param {object} root The root element to perform the search on.
 */
async function expectPeriodicTable(elements, root = screen) {
  const elementSet = new Set(elements)
  await waitFor(() => {
    elementData.elements.forEach(element => {
      const button = screen.getByText(element.symbol).closest('button')
      if (elementSet.has(element.symbol)) {
        expect(button).not.toHaveAttribute('disabled')
      } else {
        expect(button).toHaveAttribute('disabled')
      }
    })
  })
}

describe('', () => {
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/terms_aggregation_elements')
    renderSearchEntry(<InputPeriodicTable
      quantity={quantity}
      visible
    />
    )
  })
  afterEach(() => closeAPI())

  test('initial state is loaded correctly', async () => {
    // Test immediately displayed elements
    expectInputHeader(quantity)
    elementData.elements.forEach(element => {
      const name = screen.getByText(element.symbol)
      expect(name).toBeInTheDocument()
      expect(screen.getByText(element.number)).toBeInTheDocument()
      expect(screen.getByTitle(element.name)).toBeInTheDocument()
      expect(name.closest('button')).toHaveAttribute('disabled')
    })
    expect(screen.getByRole('checkbox')).toBeInTheDocument()

    // Test that only available elements are clickable after API response.
    await expectPeriodicTable(['H', 'C', 'Ti', 'Zr', 'Nb', 'Hf', 'Ta'])
  })
})

describe('', () => {
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/terms_aggregation_elements_sequence')
    renderSearchEntry(<InputPeriodicTable
      quantity={quantity}
      visible
    />
    )
  })
  afterEach(() => closeAPI())

  test('selecting an element in both non-exclusive and exclusive mode correctly updates the table', async () => {
    // Wait for hydrogen to become selectable
    await waitFor(() => { expect(screen.getByText('H').closest('button')).not.toHaveAttribute('disabled') })

    // Test that after selecting C, only H and C are selectable.
    const cButton = screen.getByText('C')
    userEvent.click(cButton)
    await expectPeriodicTable(['C', 'H'])

    // Test that after enabling exclusive search, only C is selectable
    const exclusiveCheckbox = screen.getByRole('checkbox')
    userEvent.click(exclusiveCheckbox)
    await expectPeriodicTable(['C'])
  })
})
