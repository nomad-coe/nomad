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
import {render} from "../conftest.spec"
import {QuantityItemPreview} from "./ArchiveBrowser"

test.each([
  [
    'without unit',
    3.5,
    {},
    '3.50000'],
  [
    'without display unit',
    3.5,
    {unit: 'meter'},
    '3.5·10+10 Å'
  ],
  [
    'with default display unit',
    3.5,
    {unit: 'meter', m_annotations: {display: [{unit: 'mm'}]}},
    '3500 mm'
  ]
])('Test QuantityItemPreview %s', async (name, value, def, expected) => {
  render(
    <QuantityItemPreview
      def={{name: 'value1', shape: [], type: {type_kind: 'python', type_data: 'float'}, ...def}}
      value={value}
    />
  )
  const span = document.querySelector('span')
  expect(span.textContent).toContain(expected)
})
