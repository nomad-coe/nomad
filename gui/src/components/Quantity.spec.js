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
import {render} from "./conftest.spec"
import Quantity from './Quantity'

describe('Test boolean quantity', () => {
  test.each([
    [true, 'true'],
    [false, 'false']
  ])('renders %s', (myBoolean, expectedText) => {
    const data = {
      myBoolean: myBoolean
    }
    const quantity = {
      name: 'myBoolean',
      type: { type_data: 'bool' }
    }
    const { getByText } = render(<Quantity quantity={quantity} data={data} />)
    expect(getByText(expectedText)).toBeInTheDocument()
  })
})
