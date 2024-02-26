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

describe('Test QuantityItemPreview', () => {
  it('without default display unit', async () => {
    render(
      <QuantityItemPreview
        def={{
          unit: 'meter',
          type: {type_kind: 'python', type_data: 'float'},
          shape: [],
          m_annotations: {}
        }}
        value={3.5}
      />
    )
    // should be rendered in default unit system
    const span = document.querySelector('span')
    expect(span.textContent).toContain('3.5')
    expect(span.textContent).toContain('·10')
    const sup = document.querySelector('sup')
    expect(sup.textContent).toContain('+10')
  })

  it('with default display unit', async () => {
    render(
      <QuantityItemPreview
        def={{
          unit: 'meter',
          type: {type_kind: 'python', type_data: 'float'},
          shape: [],
          m_annotations: {
            display: [{
              unit: 'mm'
            }]
          }
        }}
        value={3.5}
      />
    )
    // should be rendered in default unit provided by schema
    const span = document.querySelector('span')
    expect(span.textContent).toContain('3500 mm')
  })
})
