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
import { Quantity, parseQuantity } from './Quantity'
import { dimensionMap } from './UnitContext'

test('conversion works both ways for each compatible unit', async () => {
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
  ['same unit', 'kelvin', 'kelvin', 1, 1],
  ['temperature celsius', 'kelvin', 'celsius', 1, -272.15],
  ['temperature fahrenheit', 'kelvin', 'fahrenheit', 1, -457.87],
  ['abbreviated name', 'J', 'eV', 1, 6241509074460763000],
  ['full name', 'joule', 'electron_volt', 1, 6241509074460763000],
  ['division', 'm/s', 'angstrom/femtosecond', 1, 0.00001],
  ['multiplication', 'm*s', 'angstrom*femtosecond', 1, 9.999999999999999e+24],
  ['power with hat', 'm^2', 'angstrom^2', 1, 99999999999999980000],
  ['power with double asterisk (single)', 'm**2', 'angstrom**2', 1, 99999999999999980000],
  ['power with double asterisk (multiple)', 'm**2 / s**2', 'angstrom**2 / ms**2', 1, 99999999999999.98],
  ['explicit delta (single)', 'delta_celsius', 'delta_K', 1, 274.15],
  ['explicit delta (multiple)', 'delta_celsius / delta_celsius', 'delta_K / delta_K', 1, 1],
  ['explicit delta symbol (single)', 'Δcelsius', 'ΔK', 1, 274.15],
  ['explicit delta symbol (multiple)', 'Δcelsius / Δcelsius', 'ΔK / ΔK', 1, 1],
  ['combined', 'm*m/s^2', 'angstrom^2/femtosecond^2', 1, 9.999999999999999e-11],
  ['negative exponent', 's^-2', 'femtosecond^-2', 1, 1e-30],
  ['simple to complex with one unit', 'N', 'kg*m/s^2', 1, 1],
  ['complex to simple with one unit', 'kg*m/s^2', 'N', 1, 1],
  ['simple to complex with expression', 'N/m', 'kg/s^2', 1, 1],
  ['complex to simple with expression', 'kg/s^2', 'N/m', 1, 1],
  ['unit starting with a number', '1/minute', '1/second', 1, 0.016666666666666666]
]
)('test conversion with "to()": %s', async (name, unitA, unitB, valueA, valueB) => {
  const a = new Quantity(valueA, unitA)
  const b = a.to(unitB)
  expect(b.value()).toBeCloseTo(valueB, 10)
})

test.each([
  ['conversion with single unit', 'meter', {length: {definition: 'angstrom'}}, 1, 1e10],
  ['conversion with power', 'meter^2', {length: {definition: 'angstrom'}}, 1, 99999999999999980000],
  ['do not simplify', 'gram*angstrom/fs^2', {mass: {definition: 'kilogram'}, length: {definition: 'meter'}, time: {definition: 'second'}}, 1, 99999999999999980],
  ['do not convert to base', 'eV', {energy: {definition: 'joule'}}, 1, 1.602176634e-19],
  ['combination', 'a_u_force * angstrom', {force: {definition: 'newton'}, length: {definition: 'meter'}}, 1, 8.23872349823899e-18],
  ['use base units if derived unit not defined in system', 'newton * meter', {mass: {definition: 'kilogram'}, time: {definition: 'second'}, length: {definition: 'meter'}}, 1, 1],
  ['unit definition with prefix', 'kg^2', {mass: {definition: 'mg'}}, 1, 1e12],
  ['expression as definition', 'N', {force: {definition: '(kg m) / s^2'}}, 1, 1]
]
)('test conversion with "toSystem()": %s', async (name, unit, system, valueA, valueB) => {
  const a = new Quantity(valueA, unit)
  const b = a.toSystem(system)
  expect(b.value()).toBeCloseTo(valueB, 10)
})

test.each([
  ['celsius to kelvin', 'celsius', 'kelvin', 5, 278.15],
  ['fahrenheit to kelvin', 'fahrenheit', 'kelvin', 5, 258.15],
  ['celsius to fahrenheit', 'celsius', 'fahrenheit', 5, 41],
  ['celsius to kelvin: derived unit (implicit delta)', 'joule/celsius', 'joule/kelvin', 5, 5],
  ['celsius to kelvin: derived unit (explicit delta)', 'joule/delta_celsius', 'joule/kelvin', 5, 5],
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
  ['number only', '100', undefined, true, false, {valueString: '100', value: 100}],
  ['unit only', 'joule', null, false, true, {valueString: undefined, value: undefined, unit: new Unit('joule')}],
  ['number and unit with dimension', '100 joule', 'energy', true, true, {valueString: '100', value: 100, unit: new Unit('joule')}],
  ['number and unit without dimension', '100 joule', null, true, true, {valueString: '100', value: 100, unit: new Unit('joule')}],
  ['incorrect dimension', '100 joule', 'length', true, true, {valueString: '100', value: 100, unit: new Unit('joule'), error: 'Unit "joule" is incompatible with dimension "length"'}],
  ['missing unit', '100', 'length', true, true, {valueString: '100', value: 100, unit: undefined, error: 'Unit is required'}],
  ['missing value', 'joule', 'energy', true, true, {valueString: undefined, value: undefined, unit: new Unit('joule'), error: 'Enter a valid numerical value'}],
  ['mixing number and quantity #1', '1 / joule', 'energy^-1', false, false, {valueString: '1', value: 1, unit: new Unit('1 / joule')}],
  ['mixing number and quantity #2', '100 / joule', 'energy^-1', false, false, {valueString: '100', value: 100, unit: new Unit('1 / joule')}]

]
)('test parseQuantity: %s', async (name, input, dimension, requireValue, requireUnit, expected) => {
  const result = parseQuantity(input, dimension, requireValue, requireUnit)
  expect(result.valueString === expected.valueString).toBe(true)
  expect(result.value === expected.value).toBe(true)
  expect(result.unit?.label() === expected.unit?.label()).toBe(true)
  expect(result.error === expected.error).toBe(true)
})
