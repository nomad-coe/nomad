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

import { Quantity } from './units'
import { conversionMap } from './unitsData'

test('unit conversion works both ways for each compatible unit', async () => {
  // Create a list of all possible conversions
  const conversions = []
  for (const dimension of Object.values(conversionMap)) {
    const units = dimension.units
    for (const unitA of units) {
      for (const unitB of units) {
        conversions.push([unitA, unitB])
      }
    }
  }
  for (let [unitA, unitB] of conversions) {
    const a = new Quantity(1, unitA)
    const b = a.to(unitB)
    const c = b.to(unitA)
    expect(a.value).toBeCloseTo(c.value, 10)
  }
})

test.each([
  ['kelvin', 'kelvin'], // Same
  ['kelvin', 'celsius'], // Different
  ['J', 'eV'], // Abbreviated
  ['joule', 'electron_volt'], // Non-abbreviated
  ['m/s', 'angstrom/femtosecond'], // Division
  ['m*s', 'angstrom/femtosecond'], // Multiplication
  ['m^2', 'angstrom/femtosecond'], // Power
  ['m**2', 'angstrom/femtosecond'], // Power
  ['m*m/s^2', 'angstrom^2/femtosecond^2'] // Division+power+multiplication
]
)('complex conversions', async (unitA, unitB) => {
  const a = new Quantity(1, unitA)
  const b = a.to(unitB)
  expect(b.value).not.toBeNaN()
})

test.each([
  ['m', 'J'] // Wrong dimension
  // TODO: ['m', 'm^2'], // Incompatible expressions
]
)('invalid conversions', async (unitA, unitB) => {
  expect(() => {
    new Quantity(1, unitA).to(unitB)
  }).toThrow()
})
