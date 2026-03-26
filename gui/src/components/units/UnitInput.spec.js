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

import React, {useState} from 'react'
import { renderNoAPI, screen } from '../conftest.spec'
import userEvent from '@testing-library/user-event'
import { UnitInput } from './UnitInput'
import { unitMap } from './UnitContext'

// Find a constant and a regular unit from the test data
const constantEntry = Object.entries(unitMap).find(([, unit]) => unit.constant)
const regularEntry = Object.entries(unitMap).find(([, unit]) => !unit.constant && unit.dimension)

// UnitInput uses controlled inputValue via InputText, so we need a stateful
// wrapper that feeds typed text back as the value prop.
const StatefulUnitInput = (props) => {
  const [value, setValue] = useState('')
  return <UnitInput {...props} value={value} onChange={setValue} />
}

test('regular unit appears as a dropdown option when typed', async () => {
  renderNoAPI(<StatefulUnitInput disableGroup />)
  const input = screen.getByRole('textbox')
  await userEvent.type(input, regularEntry[0])
  const options = await screen.findAllByRole('option')
  expect(options.length).toBeGreaterThan(0)
})

test('constant does not appear as a dropdown option when typed', async () => {
  renderNoAPI(<StatefulUnitInput disableGroup />)
  const input = screen.getByRole('textbox')
  await userEvent.type(input, constantEntry[0])
  expect(screen.queryByText(constantEntry[1].label, {exact: false})).not.toBeInTheDocument()
})
