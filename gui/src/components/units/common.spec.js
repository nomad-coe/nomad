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

import { Unit } from './Unit'
import { parse } from './common'

test.each([
  ['number only', '100', {}, {valueString: '100', value: 100}],
  ['complicated number 1', '-0.015e-10', {}, {valueString: '-0.015e-10', value: -0.015e-10}],
  ['complicated number 2', '-.015E-10', {}, {valueString: '-.015E-10', value: -0.015e-10}],
  ['unit only', 'joule', {}, {valueString: undefined, value: undefined, unit: new Unit('joule')}],
  ['whitespaces', '  100  joule  ', {}, {valueString: '100', value: 100, unit: new Unit('joule') }],
  ['number and unit with dimension', '100 joule', {}, {valueString: '100', value: 100, unit: new Unit('joule')}],
  ['number and unit without dimension', '100 joule', {}, {valueString: '100', value: 100, unit: new Unit('joule')}],
  ['incorrect dimension', '100 joule', {dimension: 'length'}, {valueString: '100', value: 100, unit: new Unit('joule'), error: 'Unit "joule" is incompatible with dimension "length".'}],
  ['missing unit', '100', {requireUnit: true}, {error: 'Unit is required'}],
  ['missing unit with dimension specified', '100', {dimension: 'energy', requireUnit: true}, {error: 'Unit is required'}],
  ['missing value', 'joule', {requireValue: true}, {error: 'Enter a valid numerical value'}],
  ['mixing number and quantity #1', '1 / joule', {dimension: 'energy^-1'}, {valueString: '1', value: 1, unit: new Unit('1 / joule')}],
  ['mixing number and quantity #2', '100 / joule', {dimension: 'energy^-1'}, {valueString: '100', value: 100, unit: new Unit('1 / joule')}]

]
)('test parse: %s', async (name, input, options, expected) => {
  const result = parse(input, options)
  expect(result.valueString === expected.valueString).toBe(true)
  expect(result.value === expected.value).toBe(true)
  expect(result.unit?.label() === expected.unit?.label()).toBe(true)
  expect(result.error === expected.error).toBe(true)
})
