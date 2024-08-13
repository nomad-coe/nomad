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
import {QuantityItemPreview, QuantityValue} from "./ArchiveBrowser"

test.each([
  [
    'without unit',
    3.5,
    undefined,
    undefined,
    undefined,
    '3.50000'],
  [
    'without display unit',
    3.5,
    'meter',
    undefined,
    undefined,
    '3.5·10+10 Å'
  ],
  [
    'with default display unit',
    3.5,
    'meter',
    'mm',
    undefined,
    '3500 mm'
  ],
  [
    'complex unit with no display unit',
    3.5,
    'm**2 / second**2',
    undefined,
    undefined,
    '3.5·10-10 Å^2 / fs^2'
  ],
  [
    'complex unit with display unit',
    3.5,
    'm**2 / second**2',
    'm**2 / fs**2',
    undefined,
    '3.5·10-30 m^2 / fs^2'
  ],
  [
    'deprecated display unit in eln annotation',
    3.5,
    'm',
    undefined,
    'mm',
    '3500 mm'
  ]
])('Test QuantityItemPreview %s', async (name, value, unit, displayUnit, elnUnit, expected) => {
  const def = {
    name: 'name',
      m_def: 'nomad.metainfo.metainfo.Quantity',
      unit: unit,
      m_annotations: {
      display: displayUnit && [{
        unit: displayUnit
      }],
        eln: elnUnit && [{
        defaultDisplayUnit: elnUnit
      }]
    }
  }

  render(
    <QuantityItemPreview
      def={{name: 'value1', shape: [], type: {type_kind: 'python', type_data: 'float'}, ...def}}
      value={value}
    />
  )
  const span = document.querySelector('span')
  expect(span.textContent).toContain(expected)
})

describe("Test QuantityValue", () => {
  // eslint-disable-next-line no-undef
  const element = HTMLElement

  const originalOffsetHeight = Object.getOwnPropertyDescriptor(element.prototype, 'offsetHeight')
  const originalOffsetWidth = Object.getOwnPropertyDescriptor(element.prototype, 'offsetWidth')

  beforeAll(() => {
    Object.defineProperty(element.prototype, 'offsetHeight', { configurable: true, value: 50 })
    Object.defineProperty(element.prototype, 'offsetWidth', { configurable: true, value: 50 })
  })

  afterAll(() => {
    Object.defineProperty(element.prototype, 'offsetHeight', originalOffsetHeight)
    Object.defineProperty(element.prototype, 'offsetWidth', originalOffsetWidth)
  })

  it.each([
    [
      'without unit',
      [3.5],
      [1],
      undefined,
      undefined,
      undefined,
      '3.50000',
      false,
      '(1)',
      undefined
    ],
    [
      'without display unit',
      [3.5],
      [1],
      'meter',
      undefined,
      undefined,
      '3.5·10+10',
      true,
      '(1)',
      'Å'
    ],
    [
      'with default display unit',
      [3.5],
      [1],
      'meter',
      'mm',
      undefined,
      '3500',
      false,
      '(1)',
      'mm'
    ],
    [
      'complex unit with no display unit',
      [3.5],
      [1],
      'm**2 / second**2',
      undefined,
      undefined,
      '3.5·10-10',
      true,
      '(1)',
      'Å^2 / fs^2'
    ],
    [
      'complex unit with display unit',
      [3.5],
      [1],
      'm**2 / second**2',
      'm**2 / fs**2',
      undefined,
      '3.5·10-30',
      true,
      '(1)',
      'm^2 / fs^2'
    ],
    [
      'deprecated display unit in eln annotation',
      [3.5],
      [1],
      'm',
      undefined,
      'mm',
      '3500',
      false,
      '(1)',
      'mm'
    ]
  ])('%s', async (name, value, shape, unit, displayUnit, elnUnit, expectedValue, scientific, expectedDim, expectedUnit) => {
    const def = {
      name: 'name',
      m_def: 'nomad.metainfo.metainfo.Quantity',
      unit: unit,
      shape: shape,
      m_annotations: {
        display: displayUnit && [{
          unit: displayUnit
        }],
        eln: elnUnit && [{
          defaultDisplayUnit: elnUnit
        }]
      }
    }

    const screen = render(
      <QuantityValue
        def={{name: 'value1', type: {type_kind: 'python', type_data: 'float'}, ...def}}
        value={value}
      />
    )

    const spans = document.querySelectorAll('span')

    if (scientific) {
      const value = spans[0]
      const dim = spans[1]
      expect(value.textContent).toContain(expectedValue)
      expect(dim.textContent).toContain(expectedDim)
    } else {
      screen.getByText(expectedValue)
      const dim = spans[0]
      expect(dim.textContent).toContain(expectedDim)
    }

    if (expectedUnit) {
      screen.getByText(expectedUnit)
    }
  })
})
