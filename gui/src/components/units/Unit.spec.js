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
import { Unit } from './Unit'
import { unitMap } from './UnitContext'

test('each unit can be created using its full name, alias or short form (+ all available prefixes)', async () => {
  for (const [name, def] of Object.entries(unitMap)) {
    // Full name + prefixes
    expect(new Unit(name)).not.toBeNaN()
    if (def.prefixes) {
      for (const prefix of Object.keys(UnitMathJS.PREFIXES[def.prefixes.toUpperCase()])) {
        expect(new Unit(`${prefix}${name}`)).not.toBeNaN()
      }
    }
    // Aliases + prefixes
    if (def.aliases) {
      for (const alias of def.aliases) {
        expect(new Unit(alias)).not.toBeNaN()
        if (def.prefixes) {
          for (const prefix of Object.keys(UnitMathJS.PREFIXES[def.prefixes.toUpperCase()])) {
            expect(new Unit(`${prefix}${alias}`)).not.toBeNaN()
          }
        }
      }
    }
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
  const a = new Unit(unit)
  expect(a.label()).toBe(label)
})

test.each([
  ['dimensionless', 'dimensionless', 'dimensionless', false],
  ['single unit', 'meter', 'length', false],
  ['fixed order 1', 'meter * second', 'length time', false],
  ['fixed order 2', 'second * meter', 'length time', false],
  ['power', 'meter^3 * second^-1', 'length^3 time^-1', false],
  ['in derived', 'joule', 'energy', false],
  ['in base units', 'joule', 'mass length^2 time^-2', true]
]
)('test getting dimension": %s', async (name, unit, dimension, base) => {
  const a = new Unit(unit)
  expect(a.dimension(base)).toBe(dimension)
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
  ['power with double asterisk (single)', 'm**2', 'angstrom**2', 'Å^2'],
  ['power with double asterisk (multiple)', 'm**2 / s**2', 'angstrom**2 / ms**2', 'Å^2 / ms^2'],
  ['explicit delta (single)', 'delta_celsius', 'delta_K', 'K'],
  ['explicit delta (multiple)', 'delta_celsius / delta_celsius', 'delta_K / delta_K', 'K / K'],
  ['explicit delta symbol (single)', 'Δcelsius', 'ΔK', 'K'],
  ['explicit delta symbol (multiple)', 'Δcelsius / Δcelsius', 'ΔK / ΔK', 'K / K'],
  ['combined', 'm*m/s^2', 'angstrom^2/femtosecond^2', 'Å^2 / fs^2'],
  ['negative exponent', 's^-2', 'femtosecond^-2', 'fs^-2'],
  ['simple to complex with one unit', 'N', 'kg*m/s^2', '(kg m) / s^2'],
  ['complex to simple with one unit', 'kg*m/s^2', 'N', 'N'],
  ['simple to complex with expression', 'N/m', 'kg/s^2', 'kg / s^2'],
  ['complex to simple with expression', 'kg/s^2', 'N/m', 'N / m'],
  ['unit starting with a number', '1/minute', '1/second', 's^-1']
]
)('test conversion with "to()": %s', async (name, unitA, unitB, labelB) => {
  const a = new Unit(unitA)
  const b = a.to(unitB)
  expect(b.label()).toBe(labelB)
})

test.each([
  ['conversion with single unit', 'meter', {length: {definition: 'angstrom'}}, 'Å'],
  ['conversion with power', 'meter^2', {length: {definition: 'angstrom'}}, 'Å^2'],
  ['do not simplify', 'gram*angstrom/fs^2', {mass: {definition: 'kilogram'}, length: {definition: 'meter'}, time: {definition: 'second'}}, '(kg m) / s^2'],
  ['do not convert to base', 'eV', {energy: {definition: 'joule'}}, 'J'],
  ['combination', 'a_u_force * angstrom', {force: {definition: 'newton'}, length: {definition: 'meter'}}, 'N m'],
  ['use base units if derived unit not defined in system', 'newton * meter', {mass: {definition: 'kilogram'}, time: {definition: 'second'}, length: {definition: 'meter'}}, '(kg m m) / s^2'],
  ['unit definition with prefix', 'kg^2', {mass: {definition: 'mg'}}, 'mg^2'],
  ['expression as definition', 'N', {force: {definition: '(kg m) / s^2'}}, '(kg m) / s^2']
]
)('test conversion with "toSystem()": %s', async (name, unit, system, label) => {
  const a = new Unit(unit)
  const b = a.toSystem(system)
  expect(b.label()).toBe(label)
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
    new Unit(unitA).to(unitB)
  }).toThrow()
})
