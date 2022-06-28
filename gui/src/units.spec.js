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

import { Unit as UnitMathJS } from 'mathjs'
import { Quantity, dimensionMap, unitMap } from './units'

test('each unit can be created using its full name, alias or short form (+ all available prefixes)', async () => {
  for (const [name, def] of Object.entries(unitMap)) {
    // Full name + prefixes
    expect(new Quantity(1, name)).not.toBeNaN()
    if (def.prefixes) {
      for (const prefix of Object.keys(UnitMathJS.PREFIXES[def.prefixes.toUpperCase()])) {
        expect(new Quantity(1, `${prefix}${name}`)).not.toBeNaN()
      }
    }
    // Aliases + prefixes
    if (def.aliases) {
      for (const alias of def.aliases) {
        expect(new Quantity(1, alias)).not.toBeNaN()
        if (def.prefixes) {
          for (const prefix of Object.keys(UnitMathJS.PREFIXES[def.prefixes.toUpperCase()])) {
            expect(new Quantity(1, `${prefix}${alias}`)).not.toBeNaN()
          }
        }
      }
    }
  }
})

test('unit conversion works both ways for each compatible unit', async () => {
  // Create a list of all possible conversions
  const conversions = []
  for (const dimension of Object.values(dimensionMap)) {
    const units = dimension.units
    for (const unitA of units) {
      for (const unitB of units) {
        conversions.push([unitA, unitB])
      }
    }
  }
  for (const [unitA, unitB] of conversions) {
    const a = new Quantity(1, unitA)
    const b = a.to(unitB)
    const c = b.to(unitA)
    expect(a.value()).toBeCloseTo(c.value(), 10)
  }
})

test.each([
  ['dimensionless', 'dimensionless', ''],
  ['non-abbreviated', 'celsius', '°C'],
  ['abbreviated', '°C', '°C'],
  ['prefix long', 'millikelvin', 'mK'],
  ['prefix short', 'mK', 'mK'],
  ['division', 'angstrom / meter', 'Å / m'],
  ['multiplication', 'angstrom * meter', 'Å m'],
  ['preserve order', 'meter*second', 'm s'],
  ['preserve order', 'second*meter', 's m'],
  ['power', 'meter^2', 'm^2'],
  ['negative power', 'meter^-1', 'm^-1'],
  ['chain', 'meter*meter/second^2', '(m m) / s^2']
]
)('label abbreviation: %s', async (name, unit, label) => {
  const a = new Quantity(1, unit)
  expect(a.label()).toBe(label)
})

test.each([
  ['same unit', 'kelvin', 'kelvin', 'K'],
  ['temperature celsius', 'kelvin', 'celsius', '°C'],
  ['temperature fahrenheit', 'kelvin', 'fahrenheit', '°F'],
  ['abbreviated name', 'J', 'eV', 'eV'],
  ['full name', 'joule', 'electron_volt', 'eV'],
  ['division', 'm/s', 'angstrom/femtosecond', 'Å / fs'],
  ['multiplication', 'm*s', 'angstrom*femtosecond', 'Å fs'],
  ['power with hat', 'm^2', 'angstrom^2', 'Å^2'],
  ['power with double asterisk', 'm**2', 'angstrom**2', 'Å^2'],
  ['combined', 'm*m/s^2', 'angstrom^2/femtosecond^2', 'Å^2 / fs^2'],
  ['negative exponent', 's^-2', 'femtosecond^-2', 'fs^-2'],
  ['simple to complex with one unit', 'N', 'kg*m/s^2', '(kg m) / s^2'],
  ['complex to simple with one unit', 'kg*m/s^2', 'N', 'N'],
  ['simple to complex with expression', 'N/m', 'kg/s^2', 'kg / s^2'],
  ['complex to simple with expression', 'kg/s^2', 'N/m', 'N / m']
]
)('test conversion with "to()": %s', async (name, unitA, unitB, labelB) => {
  const a = new Quantity(1, unitA)
  const b = a.to(unitB)
  expect(b.value()).not.toBeNaN()
  expect(b.label()).toBe(labelB)
})

test.each([
  ['conversion with single unit', 'meter', {length: {name: 'angstrom'}}, 'Å'],
  ['conversion with power', 'meter^2', {length: {name: 'angstrom'}}, 'Å^2'],
  ['do not simplify', 'gram*angstrom/fs^2', {mass: {name: 'kilogram'}, length: {name: 'meter'}, time: {name: 'second'}}, '(kg m) / s^2'],
  ['do not convert to base', 'eV', {energy: {name: 'joule'}}, 'J'],
  ['combination', 'a_u_force * angstrom', {force: {name: 'newton'}, length: {name: 'meter'}}, 'N m'],
  ['use base units if derived unit not defined in system', 'newton * meter', {mass: {name: 'kilogram'}, time: {name: 'second'}, length: {name: 'meter'}}, '(m kg m) / s^2']
]
)('test conversion with "toSystem()": %s', async (name, unit, system, label) => {
  const a = new Quantity(1, unit)
  const b = a.toSystem(system)
  expect(b.value()).not.toBeNaN()
  expect(b.label()).toBe(label)
})

test.each([
  ['dimensionless', 'dimensionless', 'dimensionless'],
  ['single unit', 'meter', 'length'],
  ['fixed order 1', 'meter * second', 'length time'],
  ['fixed order 2', 'second * meter', 'length time'],
  ['power', 'meter^3 * second^-1', 'length^3 time^-1'],
  ['in derived', 'joule', 'energy', false],
  ['in base', 'joule', 'length^2 time^-2 mass']
]
)('test getting dimension": %s', async (name, unit, dimension, base = true) => {
  const a = new Quantity(1, unit)
  expect(a.dimension(base)).toBe(dimension)
})

test.each([
  ['celsius to kelvin', 'celsius', 'kelvin', 5, 278.15],
  ['fahrenheit to kelvin', 'fahrenheit', 'kelvin', 5, 258.15],
  ['celsius to fahrenheit', 'celsius', 'fahrenheit', 5, 41],
  ['celsius to kelvin: derived unit (offset not applied)', 'joule/celsius', 'joule/kelvin', 5, 5],
  ['fahrenheit to kelvin: derived unit (offset not applied)', 'joule/fahrenheit', 'joule/kelvin', 5, 9 / 5 * 5],
  ['celsius to fahrenheit: derived unit (offset not applied)', 'joule/celsius', 'joule/fahrenheit', 5, 5 / 9 * 5]
]
)('test temperature conversion": %s', async (name, unitA, unitB, valueA, valueB) => {
  const a = new Quantity(valueA, unitA)
  const b = a.to(unitB)
  const c = b.to(unitA)
  expect(b.value()).toBeCloseTo(valueB, 10)
  expect(c.value()).toBeCloseTo(valueA, 10)
})

test.each([
  [0],
  [1],
  [2],
  [3]
]
)('test different value dimensions: %sD', async (dimension) => {
  let value = 1
  for (let i = 0; i < dimension; ++i) {
    value = [value]
  }
  const a = new Quantity(value, 'angstrom')
  const b = a.to('nanometer')
  let valueA = a.value()
  let valueB = b.value()
  for (let i = 0; i < dimension; ++i) {
    valueA = valueA[0]
    valueB = valueB[0]
  }
  expect(valueA).toBeCloseTo(10 * valueB)
})

test.each([
  ['incompatible dimensions', 'm', 'J'],
  ['wrong power of the correct dimension', 'm', 'm^2'],
  ['unknown unit', 'm', 'beard-second'],
  ['unknown operator', 'm%2', 'Å%2'],
  ['unclosed bracket', '(m^2', '(Å^2']
]
)('invalid conversions with "to()": %s', async (name, unitA, unitB) => {
  expect(() => {
    new Quantity(1, unitA).to(unitB)
  }).toThrow()
})
