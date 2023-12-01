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

import {isNumber, isArray, isNil} from 'lodash'
import {Unit} from './Unit'
import {mapDeep} from '../../utils'

/**
 * Class for persisting persisting a numeric value together with unit
 * information.
 */
export class Quantity {
  /**
   * @param {number | n-dimensional array of numbers} value Numeric value. See
   * also the argument 'normalized'.
   * @param {str | Unit} unit Unit for the quantity.
   * @param {boolean} normalized Whether the given numeric value is already
   * normalized to base units.
   */
  constructor(value, unit, normalized = false) {
    this.unit = new Unit(unit)
    if (!isNumber(value) && !isArray(value)) {
      throw Error('Please provide the the value as a number, or as a multidimensional array of numbers.')
    }

    // This attribute stores the quantity value in 'normalized' form that is
    // given in the base units (=SI). This value should only be determined once
    // during the unit initialization and all calls to value() will then lazily
    // determine the value in the currently set units. This avoids 'drift' in
    // the value caused by several consecutive changes of the units.
    this.normalized_value = normalized ? value : this.normalize(value)
  }

  /**
   * Get value in current units.
   * @returns The numeric value in the currently set units.
   */
  value() {
    return this.denormalize(this.normalized_value)
  }

  /**
   * Convert value from currently set units to base units.
   * @param {n-dimensional array} value Value in currently set units.
   * @returns Value in base units.
   */
  normalize(value) {
    return mapDeep(value, (x) => this.unit.mathjsUnit._normalize(x))
  }

  /**
   * Convert value from base units to currently set units.
   * @param {n-dimensional array} value Value in base units.
   * @returns Value in currently set units.
   */
  denormalize(value) {
    return mapDeep(value, (x) => this.unit.mathjsUnit._denormalize(x))
  }

  label() {
    return this.unit.label()
  }

  dimension(base) {
    return this.unit.dimension(base)
  }

  to(unit) {
    return new Quantity(this.normalized_value, this.unit.to(unit), true)
  }

  toSI() {
    return new Quantity(this.normalized_value, this.unit.toSI(), true)
  }

  toSystem(system) {
    return new Quantity(this.normalized_value, this.unit.toSystem(system), true)
  }

  /**
   * Checks if the given Quantity is equal to this one.
   * @param {Quantity} quantity Quantity to compare to
   * @returns boolean Whether quantities are equal
   */
  equal(quantity) {
    if (quantity instanceof Quantity) {
      return this.normalized_value === quantity.normalized_value && this.unit.equalBase(quantity.unit)
    } else {
      throw Error('The given value is not an instance of Quantity.')
    }
  }
}

/**
 * Convenience function for parsing value and unit information from a string.
 *
 * @param {string} input The input string to parse
 * @param {boolean} requireValue Whether a value is required.
 * @param {boolean} requireUnit Whether a unit is required.
 * @param {string} dimension Dimension for the unit. Nil value means a
 * dimensionless unit.
 * @returns Object containing the following properties, if available:
 *  - value: Numerical value as a number
 *  - valueString: Numerical value as a string
 *  - unit: Unit instance
 *  - unitString: Unit as a string
 *  - error: Error messsage
 */
export function parseQuantity(input, requireValue = true, requireUnit = true, dimension = undefined) {
  input = input.trim()
  const valueString = input.match(/^[+-]?((\d+\.\d+|\d+\.|\.\d?|\d+)(e|e\+|e-)\d+|(\d+\.\d+|\d+\.|\.\d?|\d+))?/)?.[0]
  if (requireValue && isNil(valueString)) {
    return {error: 'Enter a valid numerical value'}
  }
  const value = Number(valueString)
  const unitString = input.substring(valueString.length).trim()
  const dim = isNil(dimension) ? 'dimensionless' : dimension
  if (unitString === '' && dim !== 'dimensionless' && requireUnit) {
    return {value, valueString, unitString, error: 'Unit is required'}
  }
  if (unitString === '' && !requireUnit) {
    return {value, valueString, unitString}
  }
  if (dim === 'dimensionless' && unitString !== '') {
    return {value, valueString, unitString, error: 'Enter a numerical value without units'}
  }
  let unit
  try {
    unit = new Unit(dim === 'dimensionless' ? 'dimensionless' : input)
  } catch {
    return {valueString, value, unitString, error: `Unit "${unitString}" is not available`}
  }
  const inputDim = unit.dimension(false)
  if (inputDim !== dimension) {
    return {valueString, value, unitString, unit, error: `Unit "${unitString}" has incompatible dimension`}
  }
  return {value, valueString, unit, unitString}
}
