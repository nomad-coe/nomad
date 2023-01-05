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
import {
  closeAPI,
  render,
  screen, startAPI, wait, waitForGUI
} from '../conftest.spec'
import {EditQuantityExamples} from './EditQuantityExamples'
import {within} from '@testing-library/dom'
import {fireEvent, waitFor} from '@testing-library/react'

const changeValue = (input, newValue) => {
  fireEvent.change(input, { target: { value: newValue } })
  fireEvent.keyDown(input, {key: 'Enter', code: 'Enter'})
}

const expectValue = async (input, expectedValue, expectedStoredValue = undefined) => {
  await waitFor(() => expect(input.value).toEqual(expectedValue))
  if (expectedStoredValue) {
    await waitFor(() => expect(screen.queryByText(expectedStoredValue)).toBeInTheDocument())
  }
}

test('correctly renders edit quantities', async () => {
  render(<EditQuantityExamples />)

  // Wait to load the entry metadata, i.e. wait for some texts to appear
  await screen.findByText('three')

  const numberFieldValue = screen.queryAllByTestId('number-edit-quantity-value')
  const numberFieldValueInput = within(numberFieldValue[1]).getByRole('textbox')
  const numberFieldUnit = screen.queryAllByTestId('number-edit-quantity-unit')
  const numberFieldUnitInput = within(numberFieldUnit[0]).getByRole('textbox', { hidden: true })
  const numberFieldCheckBoxs = screen.queryAllByTitle('If checked, numeric value is converted when the unit is changed.')
  const numberFieldCheckBox = within(numberFieldCheckBoxs[0]).getByRole('checkbox', { hidden: true })
  await waitFor(() => expect(numberFieldUnitInput.value).toEqual('fs'))

  // Check that 'read' mode is enabled
  await waitFor(() => expect(numberFieldCheckBox.checked).toBe(true))

  // Write value in text field, press enter, see that normalize value is written
  // in debug output and that the value in the text field does not change (the
  // value that is returned has gone through some conversions and due to
  // floating point accuracies it's serialized form may change, but detect this and prevent
  // the text from changing)
  await changeValue(numberFieldValueInput, '1')
  await expectValue(numberFieldValueInput, '1', /"float_unit": 1e-15/i)

  // Change the unit, see that text input changes, debug output remains the same
  await changeValue(numberFieldUnitInput, 's')
  await expectValue(numberFieldValueInput, '1e-15', /"float_unit": 1e-15/i)

  // Enter value with unit, see that only numeric value is preserved in field,
  // unit selection has changed and that debug output is correct
  await changeValue(numberFieldValueInput, '30000000fs')
  await expectValue(numberFieldValueInput, '30000000', /"float_unit": 3e-8/i)
  await expectValue(numberFieldUnitInput, 'fs')

  await changeValue(numberFieldValueInput, '1minute')
  await expectValue(numberFieldValueInput, '1', /"float_unit": 60/i)
  await expectValue(numberFieldUnitInput, 'minute')

  // Change the unit, see that text input changes, debug output remains the same
  await changeValue(numberFieldUnitInput, 'fs')
  await expectValue(numberFieldValueInput, '6e+16', /"float_unit": 60/i)

  // Change to 'write' mode
  fireEvent.click(numberFieldCheckBox)
  await wait(undefined, 100)
  await waitFor(() => expect(numberFieldCheckBox.checked).toBe(false))

  // Change the unit, see that text input remains the same, debug output changes
  await changeValue(numberFieldUnitInput, 's')
  await expectValue(numberFieldValueInput, '6e+16', /"float_unit": 60000000000000000/i)

  const numberFieldValueInputInMeter = within(numberFieldValue[4]).getByRole('textbox')
  const numberFieldUnitInputInMeter = within(numberFieldUnit[3]).getByRole('textbox', { hidden: true })

  // Enter value with unit, see that only numeric value is preserved in field,
  // unit selection has changed and that debug output is correct
  await changeValue(numberFieldValueInputInMeter, '10.5angstrom')
  await expectValue(numberFieldValueInputInMeter, '10.5', /"float_with_bounds": 1\.05e-9/i)
  await expectValue(numberFieldUnitInputInMeter, 'Å')

  // check maximum minimum bounds during unit changes
  await changeValue(numberFieldValueInputInMeter, '10.5m')
  await waitFor(() => expect(screen.queryByText(/Enter a value that is equal or smaller than 10 \(meter\)/i)).toBeInTheDocument())

  await changeValue(numberFieldValueInputInMeter, '10m')
  await expectValue(numberFieldValueInputInMeter, '10', /"float_with_bounds": 10/i)
  await expectValue(numberFieldUnitInputInMeter, 'm')

  // check maximum minimum bounds without changing unit
  await changeValue(numberFieldValueInputInMeter, '-1')
  await waitFor(() => expect(screen.queryByText(/Enter a value that is equal or larger than 0 \(meter\)/i)).toBeInTheDocument())

  // Test for the URLEditQuantity
  const UrlComponent = screen.getByTestId('URLEditQuantity')
  const invalidUrlMsg = () => within(UrlComponent).queryByText(/invalid url string!/i)
  const UrlTextbox = () => within(UrlComponent).getByRole('textbox')
  const redirectButton = () => within(UrlComponent).queryByRole('button', { name: /open_link/i })

  expect(invalidUrlMsg()).not.toBeInTheDocument()
  expect(redirectButton()).not.toBeInTheDocument()
  fireEvent.change(UrlTextbox(), { target: { value: 'a' } })
  await waitFor(() => expect(invalidUrlMsg()).toBeInTheDocument())

  fireEvent.change(UrlTextbox(), { target: { value: 'https://nomad-lab.eu/' } })
  await waitFor(() => expect(invalidUrlMsg()).not.toBeInTheDocument())
  await waitFor(() => expect(redirectButton()).toBeInTheDocument())
})

test('Test AuthorEditQuantity', async () => {
  await startAPI('tests.states.uploads.empty', 'tests/data/editquantity/user')
  render(<EditQuantityExamples />)

  // Wait to load the entry metadata, i.e. wait for some texts to appear
  await screen.findByText('User account')

  const userFields = screen.getAllByTestId('user-edit-quantity')
  expect(userFields.length).toBe(2)

  const userField = userFields[1]
  const userFieldInput = within(userField).getByRole('textbox')
  userField.focus()
  // assign an incomplete value to the input field
  fireEvent.change(userFieldInput, { target: { value: 'schei' } })
  await waitForGUI(700, true)
  await waitFor(() => expect(userFieldInput.value).toEqual('schei'))
  await waitForGUI(700, true)
  fireEvent.keyDown(userField, { key: 'ArrowDown' })
  fireEvent.keyDown(userField, { key: 'Enter' })
  await waitForGUI()
  await waitFor(() => expect(userFieldInput.value).toEqual('Markus Scheidgen (FHI)'))

  const firstNameField = screen.getByText('First name').parentElement
  const firstNameFieldInput = within(firstNameField).getByRole('textbox')
  await waitFor(() => expect(firstNameFieldInput.value).toEqual('Markus'))

  const lastNameField = screen.getByText('Last name').parentElement
  const lastNameFieldInput = within(lastNameField).getByRole('textbox')
  await waitFor(() => expect(lastNameFieldInput.value).toEqual('Scheidgen'))

  const affiliationField = screen.getByText('Affiliation').parentElement
  const affiliationFieldInput = within(affiliationField).getByRole('textbox')
  await waitFor(() => expect(affiliationFieldInput.value).toEqual('FHI'))

  const emailField = screen.getByText('Email').parentElement
  const emailFieldInput = within(emailField).getByRole('textbox')
  await waitFor(() => expect(emailFieldInput.value).toEqual(''))

  const addressField = screen.getByText('Address').parentElement
  const addressFieldInput = within(addressField).getByRole('textbox')
  await waitFor(() => expect(addressFieldInput.value).toEqual('Berlin'))

  closeAPI()
})
