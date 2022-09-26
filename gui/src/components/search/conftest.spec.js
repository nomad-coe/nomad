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
import PropTypes from 'prop-types'
import assert from 'assert'
import { waitFor } from '@testing-library/dom'
import elementData from '../../elementData.json'
import { screen, WrapperDefault, WrapperNoAPI } from '../conftest.spec'
import { render } from '@testing-library/react'
import { SearchContext } from './SearchContext'
import { filterData } from './FilterRegistry'
import { format } from 'date-fns'
import { DType } from '../../utils'

/*****************************************************************************/
// Renders
/**
 * Entry search render.
 */
const WrapperSearch = ({children}) => {
  return <WrapperDefault>
    <SearchContext resource="entries">
      {children}
    </SearchContext>
  </WrapperDefault>
}

WrapperSearch.propTypes = {
  children: PropTypes.node
}

export const renderSearchEntry = (ui, options) =>
  render(ui, {wrapper: WrapperSearch, ...options})

/**
 * Entry search render without API.
 */
const WrapperNoAPISearch = ({children}) => {
  return <WrapperNoAPI>
    <SearchContext resource="entries">
      {children}
    </SearchContext>
  </WrapperNoAPI>
}

WrapperNoAPISearch.propTypes = {
  children: PropTypes.node
}

export const renderNoAPISearchEntry = (ui, options) =>
  render(ui, {wrapper: WrapperNoAPISearch, ...options})

/*****************************************************************************/
// Expects
/**
 * Tests that the initial state of an InputHeader is correct.
 *
 * @param {string} Quantity Full metainfo name for the quantity.
 * @param {boolean} disableScale Is the statistics scaling is disabled.
 * @param {object} root The container to work on.
 */
export async function expectInputHeader(quantity, disableScale, root = screen) {
  const data = filterData[quantity]
  const label = data.label
  const description = data.description
  expect(await root.findByText(label, {exact: false})).toBeInTheDocument()
  expect(root.getByTitle(description)).toBeInTheDocument()
  if (!disableScale) {
    const scale = data.scale
    expect(root.getByText(scale)).toBeInTheDocument()
  }
}

/**
 * Tests that an InputList is rendered with the given contents.
 * @param {string} quantity The quantity name
 * @param {bool} loaded Whether the data is already loaded.
 * @param {string[]} items List of items to be displayed
 * @param {string} prompt The prompt to show at the end. One of 'all', 'first'.
 * If the given list of items is empty, this prompt is ignored.
 */
export async function expectInputList(quantity, loaded, items, prompt, root = screen) {
    const prompts = new Set(['all', 'top', 'single'])
    assert(
      items.length === 0 || prompts.has(prompt),
      `Please provide one of the values: ${[...prompts].join(', ')}`
    )
    assert(
      prompt !== 'single' || items.length === 1,
      'Only provide one value with prompt=single'
    )

    // Test immediately displayed elements
    await expectInputHeader(quantity)

    // Check that placeholder disappears
    if (!loaded) {
      await waitFor(() => expect(root.queryByTestId('inputlist-placeholder')).toBe(null))
    }

    // Test elements that are displayed after API response
    for (const item of items) {
      expect(await root.findByText(item)).toBeInTheDocument()
    }

    // Expect a message at the end
    if (items.length === 0) {
      expect(root.getByText('No options available with current query.')).toBeInTheDocument()
    } else if (prompt === 'all') {
      expect(root.getByText(`Showing all ${items.length} items`)).toBeInTheDocument()
    } else if (prompt === 'top') {
      expect(root.getByText(`Showing top ${items.length} items`)).toBeInTheDocument()
    } else if (prompt === 'single') {
      expect(root.getByText(`Showing the only item`)).toBeInTheDocument()
    }
}

/**
 * Tests that an InputRange is rendered with the given contents.
 * @param {string} quantity The quantity name
 * @param {bool} loaded Whether the data is already loaded
 * @param {bool} histogram Whether the histogram is shown
 * @param {bool} placeholder Whether the placeholder should be checked
 */
export async function expectInputRange(quantity, loaded, histogram, anchored, min, max, root = screen) {
    // Test header
    await expectInputHeader(quantity, true)

    // Check histogram
    if (histogram) {
      // Check that placeholder disappears
      if (!loaded) {
        await waitFor(() => expect(root.queryByTestId('inputrange-histogram-placeholder')).toBe(null))
      }
    }

    // Test text elements if the component is not anchored
    if (!anchored) {
      const data = filterData[quantity]
      const dtype = data.dtype
      if (dtype === DType.Timestamp) {
        expect(root.getByText('Start time')).toBeInTheDocument()
        expect(root.getByText('End time')).toBeInTheDocument()
      } else {
        expect(root.getByText('min')).toBeInTheDocument()
        expect(root.getByText('max')).toBeInTheDocument()
      }

      // Get the formatted datetime in current timezone (timezones differ, so the
      // local timezone must be used in order to prevent tests from breaking).
      if (dtype === DType.Timestamp) {
        min = format(min, 'dd/MM/yyyy kk:mm')
        max = format(max, 'dd/MM/yyyy kk:mm')
      }

      // Test elements that are displayed after API response
      expect(await root.findByDisplayValue(min)).toBeInTheDocument()
      expect(await root.findByDisplayValue(max)).toBeInTheDocument()
    }
}

/**
 * Tests that an InputPeriodicTable is rendered with the given contents.
 * @param {string} quantity The quantity name
 * @param {bool} loaded Whether the data is already loaded.
 * @param {array} elements List of chemical symbols.
 * @param {object} root The root element to perform the search on.
 */
export async function expectInputPeriodicTable(quantity, loaded, elements, root = screen) {
    // Test that all elements are displayed
    elementData.elements.forEach(element => {
      const name = root.getByText(element.symbol)
      expect(name).toBeInTheDocument()
      expect(root.getAllByText(element.number)) // This number may also be used as a count
      expect(root.getByTitle(element.name)).toBeInTheDocument()
      if (!loaded) {
        expect(name.closest('button')).toHaveAttribute('disabled')
      }
    })
    expect(screen.getByRole('checkbox')).toBeInTheDocument()

    // Test header
    await expectInputHeader(quantity)

    // Test that only available elements are clickable after API response.
    await expectInputPeriodicTableItems(elements)
}

/**
 * Tests that an InputPeriodicTable has the given elements available.
 * @param {array} elements List of chemical symbols.
 * @param {object} root The root element to perform the search on.
 */
export async function expectInputPeriodicTableItems(elements, root = screen) {
    // Test that only available elements are clickable after API response.
    const elementSet = new Set(elements)
    await waitFor(() => {
      elementData.elements.forEach(element => {
        const button = root.getByText(element.symbol).closest('button')
        expect(button).not.toBe(null)
        if (elementSet.has(element.symbol)) {
          expect(button).not.toHaveAttribute('disabled')
        } else {
          expect(button).toHaveAttribute('disabled')
        }
      })
    })
}
