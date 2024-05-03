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
import { within, waitFor } from '@testing-library/dom'
import elementData from '../../elementData.json'
import { screen, WrapperDefault } from '../conftest.spec'
import { render } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { SearchContext } from './SearchContext'
import { defaultFilterData } from './FilterRegistry'
import { format } from 'date-fns'
import { DType, getDisplayLabel, parseJMESPath } from '../../utils'
import { Unit } from '../units/Unit'
import { ui } from '../../config'
import { menuMap } from './menus/FilterMainMenu'

/*****************************************************************************/
// Renders
/**
 * Render within a search context.
 */
export const WrapperSearch = ({children}) => {
  return <WrapperDefault>
    <SearchContext resource="entries" id='entries'>
      {children}
    </SearchContext>
  </WrapperDefault>
}

WrapperSearch.propTypes = {
  children: PropTypes.node
}

export const renderSearchEntry = (ui, options) =>
  render(ui, {wrapper: WrapperSearch, ...options})

/*****************************************************************************/
// Expects

/**
 * Tests that the initial state of a FilterTitle is correct.
 *
 * @param {string} quantity Full metainfo name for the quantity.
 * @param {string} label Expected label.
 * @param {string} description Expected description.
 * @param {string} unit Expected unit label.
 * @param {string} disableUnit Whether to disable showing unit.
 * @param {object} root The container to work on.
 */
export async function expectFilterTitle(quantity, label, description, unit, disableUnit, root = screen) {
  const data = defaultFilterData[quantity]
  let finalLabel = label || data?.label
  const finalDescription = description || data?.description
  if (!disableUnit) {
    const finalUnit = unit || (
      data?.unit && new Unit(data?.unit).toSystem(ui.unit_systems.options.Custom.units).label()
    )
    if (finalUnit) finalLabel = `${finalLabel} (${finalUnit})`
  }
  await root.findAllByText(finalLabel)
  expect(root.getAllByTooltip(finalDescription)[0]).toBeInTheDocument()
}

/**
 * Tests that the initial state of an InputHeader is correct.
 *
 * @param {string} Quantity Full metainfo name for the quantity.
 * @param {boolean} disableScale Is the statistics scaling is disabled.
 * @param {object} root The container to work on.
 */
export async function expectInputHeader(quantity, disableScale, root = screen) {
  await expectFilterTitle(quantity)
  const data = defaultFilterData[quantity]
  if (!disableScale) {
    const scale = data.scale
    expect(root.getByText(scale)).toBeInTheDocument()
  }
}

/**
 * Tests that a WidgetTerms is rendered with the given contents.
 * @param {object} widget The widget setup
 * @param {bool} loaded Whether the data is already loaded
 * @param {string[]} items List of items to be displayed
 * @param {string} prompt The prompt to show at the end. One of 'all', 'first'.
 * If the given list of items is empty, this prompt is ignored.
 */
export async function expectWidgetTerms(widget, loaded, items, prompt, root = screen) {
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
    await expectFilterTitle(widget.quantity)

    // Check that placeholder disappears
    if (!loaded) {
      await waitFor(() => expect(root.queryByTestId('widgetterms-placeholder')).toBe(null))
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
 * Tests that a WidgetScatterPlot is rendered with the given contents.
 * @param {object} widget The widget setup
 * @param {bool} loaded Whether the data is already loaded
 */
export async function expectWidgetScatterPlot(widget, loaded, colorTitle, legend, root = screen) {
    // Test immediately displayed elements
    const {quantity: x} = parseJMESPath(widget.x.quantity)
    const {quantity: y} = parseJMESPath(widget.y.quantity)
    await expectFilterTitle(x)
    await expectFilterTitle(y)
    if (colorTitle) {
      if (colorTitle.quantity) {
        await expectFilterTitle(colorTitle.quantity)
      } else {
        await expectFilterTitle(undefined, colorTitle.title, colorTitle.unit)
      }
    }
    for (const label of legend || []) {
      root.getByText(label)
    }

    // Check that placeholder disappears
    if (!loaded) {
      await waitFor(() => expect(root.queryByTestId(`${widget.id}-placeholder`)).toBe(null))
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
      const data = defaultFilterData[quantity]
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
export async function expectPeriodicTable(quantity, loaded, elements, root = screen) {
    // Test that all elements are displayed
    elementData.elements.forEach(element => {
      const name = root.getByText(element.symbol)
      expect(name).toBeInTheDocument()
      expect(root.getAllByText(element.number)) // This number may also be used as a count
      expect(root.getByTitle(element.name)).toBeInTheDocument()
      if (!loaded) {
        expectElement(element.name, true)
      }
    })
    expect(screen.getByRole('checkbox')).toBeInTheDocument()

    // Test header
    await expectInputHeader(quantity)

    // Test that only available elements are clickable after API response.
    await expectPeriodicTableItems(elements)
}

/**
 * Tests that a PeriodicTable element is displayed correctly.
 * @param {string} name Full name of the element.
 * @param {bool} disabled Whether the element should be disabled.
 */
export function expectElement(name, disabled) {
    const rect = screen.getByTestId(name)
    expect(rect).not.toBe(null)
    const classes = [...rect.classList].join(' ')
    if (disabled) {
      expect(classes).toMatch(/rectDisabled/)
    } else {
      expect(classes).not.toMatch(/rectDisabled/)
    }
}

/**
 * Tests that an InputPeriodicTable has the given elements available.
 * @param {array} elements List of chemical symbols.
 * @param {object} root The root element to perform the search on.
 */
export async function expectPeriodicTableItems(elements, root = screen) {
    // Because the component does not know when it is loading data, we simply
    // have to wait if it eventually fullfills the expected results.
    const elementSet = new Set(elements)
    await waitFor(() => {
      elementData.elements.forEach(element => {
        expectElement(element.name, !elementSet.has(element.symbol))
      })
    })
}

/**
 * Tests that the correct FilterMainMenu items are displayed.
 * @param {object} context The used search context
 * @param {object} root The root element to perform the search on.
 */
export async function expectFilterMainMenu(context, root = screen) {
    // Check that menu item labels are displayed
    const menuConfig = context.filter_menus
    const menuItems = (menuConfig.include || Object.keys(menuConfig.options))
      .filter(key => !menuConfig?.exclude?.includes(key))
      .map(key => ({key, ...menuConfig.options[key]}))
    for (const menuItem of menuItems) {
      const label = menuItem.label
      if (label) {
        const labelElement = screen.getByTestId(`menu-item-label-${menuItem.key}`)
        expect(labelElement).toBeInTheDocument()
        expect(within(labelElement).getByText(label)).toBeInTheDocument()
      }
    }

    // Check that action labels are displayed
    const actionItems = []
    for (const menuItem of menuItems) {
      for (const key of menuItem?.actions?.include || []) {
        actionItems.push(menuItem.actions.options[key])
      }
    }
    for (const actionItem of actionItems) {
      const actionLabel = actionItem.label
      if (actionLabel) {
        const labelElement = screen.getByText(actionLabel)
        expect(labelElement).toBeInTheDocument()
        expect(within(labelElement).getByText(actionLabel)).toBeInTheDocument()
      }
    }

    // Check that clicking the menu items with a submenu opens up the menu
    for (const menuItem of menuItems) {
      if (menuMap[menuItem.key]) {
        const labelMenu = screen.getByTestId(`menu-item-label-${menuItem.key}`)
        const labelSubMenu = await screen.findByTestId(`filter-menu-header-${menuItem.key}`)
        expect(labelSubMenu).not.toBeVisible()
        await userEvent.click(labelMenu)
        expect(labelSubMenu).toBeVisible()
      }
    }
}

/**
 * Tests that the correct SearchResults are displayed.
 * @param {object} context The used search context
 * @param {object} root The root element to perform the search on.
 */
export async function expectSearchResults(context, root = screen) {
    // Wait until search results are in
    expect(await screen.findByText("search result", {exact: false})).toBeInTheDocument()

    // We need to focus the search on the table itself because the column labels
    // are often found on the menu as well, and a single selector is not enough
    // to specify the target
    const container = within(screen.getByTestId('search-results'))

    // Check that correct columns are displayed
    const columnConfig = context.columns
    const columnLabels = columnConfig.selected.map(key => {
      const config = columnConfig.options[key]
      const unit = config.unit
      const label = config.label || defaultFilterData[key]?.label || getDisplayLabel({name: key.split('.').slice(-1)[0]})
      return unit
        ? `${label} (${new Unit(unit).label()})`
        : label
    })
    for (const columnLabel of columnLabels) {
      expect(container.getByText(columnLabel)).toBeInTheDocument()
    }
}
