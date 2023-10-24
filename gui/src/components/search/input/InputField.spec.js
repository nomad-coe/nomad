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
import { waitFor, within, waitForElementToBeRemoved } from '@testing-library/dom'
import { startAPI, closeAPI, screen } from '../../conftest.spec'
import { renderSearchEntry, expectInputHeader } from '../conftest.spec'
import { defaultFilterData } from '../FilterRegistry'
import InputField from './InputField'
import userEvent from '@testing-library/user-event'

const stateName = 'tests.states.search.search'
const optionsStructural = ['bulk', '2D', 'molecule / cluster']
const optionsProgramName = ['VASP', 'exciting', 'ABACUS', 'ABINIT', 'AFLOW']
const optionsXC = [
  'GGA_C_PBE_SOL',
  'GGA_X_PBE_SOL',
  'LDA_C_PZ',
  'LDA_X_PZ'
]

describe('', () => {
  const quantity = 'results.material.structural_type'
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/inputfield-structural-type')
    renderSearchEntry(<InputField
      quantity={quantity}
      disableSearch
      visible
    />)
  })
  afterEach(() => closeAPI())

  test('initial state is loaded correctly for quantity with fixed options', async () => {
    // Test immediately displayed elements
    const allOptions = getAllOptions(quantity)
    await expectInputHeader(quantity)
    for (const option of allOptions) {
      expect(await screen.findByText(option)).toBeInTheDocument()
    }

    // Test that options become selectable after API call finishes
    await expectOptions(optionsStructural, allOptions)
  })
})

describe('', () => {
  const quantity = 'results.method.simulation.program_name'
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/inputfield-program-name')
    renderSearchEntry(<InputField
      quantity={quantity}
      visible
      data-testid="inputfield"
    />)
  })
  afterEach(() => closeAPI())

  test('initial state is loaded correctly for quantity with dynamically loaded options', async () => {
    // Test that placeholder is shown while loading
    const placeholder = screen.queryByTestId('inputfield-placeholder')
    expect(placeholder).toBeInTheDocument()

    // Check that placeholder disappears
    await waitForElementToBeRemoved(() => screen.queryByTestId('inputfield-placeholder'))

    // Test header
    await expectInputHeader(quantity)

    // Test that options become selectable after API call finishes
    await expectOptions(['VASP', 'exciting'], optionsProgramName)

    // Test that the "show more" button is shown, but "show less" is not shown
    expect(screen.getByText('Show more')).toBeInTheDocument()
    expect(screen.queryByText('Show less')).not.toBeInTheDocument()
  })
})

describe('', () => {
  const quantity = 'results.material.structural_type'
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/inputfield-structural-type-edit')
    renderSearchEntry(<InputField
      quantity={quantity}
      visible
      data-testid="inputfield"
    />)
  })
  afterEach(() => closeAPI())

  test('selecting an option for an exlusive filter does not update the displayed options', async () => {
    // Test that options are selectable after API call finishes
    await expectOptions(optionsStructural, optionsStructural)

    // Select 'bulk' and test that all options are still available.
    const checkbox = queryByInputItemName('bulk')
    await userEvent.click(checkbox)
    await expectOptions(optionsStructural, optionsStructural)
  })
})

describe('', () => {
  const quantity = 'results.method.simulation.dft.xc_functional_names'
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/inputfield-xc-functional-names')
    renderSearchEntry(<InputField
      quantity={quantity}
      visible
      data-testid="inputfield"
    />)
  })
  afterEach(() => closeAPI())

  test('selecting an option for a non-exlusive filter updates the displayed options', async () => {
    // Test that options are selectable after API call finishes
    await expectOptions(optionsXC, optionsXC)

    // Select PBE exchange and test that only PBE correlation is shown
    // afterwards
    const checkbox = queryByInputItemName('GGA_C_PBE_SOL')
    await userEvent.click(checkbox)
    await expectOptions(['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'], optionsXC, false)
  })
})

describe('', () => {
  const quantity = 'results.method.simulation.dft.xc_functional_names'
  beforeEach(async () => {
    await startAPI(stateName, 'tests/data/search/inputfield-xc-functional-names-suggestion')
    renderSearchEntry(<InputField
      quantity={quantity}
      visible
      data-testid="inputfield"
    />)
  })
  afterEach(() => closeAPI())

  test('search field is shown, suggestions are shown when typing, selecting suggested values works correctly', async () => {
    // Wait for initial render
    await expectOptions(optionsXC, optionsXC, false)

    // See that the input field is shown and is wait until it is not disabled
    const input = screen.getByPlaceholderText('Type here')
    expect(input).toBeInTheDocument()
    await waitFor(() => expect(input).not.toBeDisabled())

    // Start typing a value, it is expected that a list of suggestions will be
    // shown shortly afterwards
    await userEvent.type(input, 'SOL')
    const suggestions = ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL']
    for (const suggestion of suggestions) {
      await screen.findByMenuItem(suggestion)
    }

    // Select a suggested option by clicking it. This should update the list.
    const suggestion = screen.getByMenuItem('GGA_C_PBE_SOL')
    await userEvent.click(suggestion)
    await expectOptions(suggestions, suggestions, false)
  })
})

/**
 * Tests that only the given options are selectable within an InputField.
 * @param {array} selectable Selectable options
 * @param {array} all All options
 * @param {boolean} visible Whether the non-available options should still be
 * displayed in a disabled state.
 * @param {object} root The root element to perform the search on.
 */
async function expectOptions(selectable, all, visible = true, root = screen) {
  const availableSet = new Set(selectable)
  await waitFor(() => {
    for (const option of all) {
      const inputCheckbox = queryByInputItemName(option)
      if (availableSet.has(option)) {
        expect(inputCheckbox).not.toHaveAttribute('disabled')
      } else if (visible) {
        expect(inputCheckbox).toHaveAttribute('disabled')
      } else {
        expect(inputCheckbox).toBe(null)
      }
    }
  })
}

/**
 * Finds the checkbox corresponding to an InputItem with the given value.
 * @param {string} name The option value that is displayed
 * @returns {element} The checkbox input HTML element.
 */
function queryByInputItemName(option, root = screen) {
  const inputLabel = root.queryByText(option)
  const inputCheckbox = inputLabel && within(inputLabel.closest('label')).getByRole('checkbox')
  return inputCheckbox
}

/**
 * Returns the list of all the enumerated options that are available for a
 * quantity.
 * @param {string} quantity Quantity name
 * @returns {array} List of options for the given quantity.
 */
function getAllOptions(quantity) {
  return [...Object.values(defaultFilterData[quantity].options)].map(option => option.label)
}
