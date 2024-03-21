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

import {isNumber, isArray} from 'lodash'
import {Unit, normalizeExpression} from './Unit'
import { Unit as UnitMathJS } from 'mathjs'
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
      throw Error('Please provide the value as a number, or as a multidimensional array of numbers.')
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

  dimension(base = false) {
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
 * @param {string} dimension Dimension for the unit. Note that you should use
 *   the base dimensions which you can get e.g. with .dimension(true). Defaults
 *   to 'dimensionless' if not specified. If you want to disable dimension
 *   checks, use null.
 * @param {boolean} requireValue Whether an explicit numeric value is required at the start of the input.
 * @param {boolean} requireUnit Whether an explicit unit in the input is required at the end of the input.
 * @returns Object containing the following properties, if available:
 *  - value: Numerical value as a number
 *  - valueString: The original number input as a string. Note that this can only return
 *    the number when it is used as a prefix, and does not work with numbers that are
 *    part of a complex expression, e.g. 300 eV / 1000 K.
 *  - unit: Unit instance
 *  - error: Error messsage
 */
export function parseQuantity(input, dimension = 'dimensionless', requireValue = false, requireUnit = false) {
  input = input.trim()
  let error
  let value
  let valueString = input.match(/^[+-]?((\d+\.\d+|\d+\.|\.\d?|\d+)(e|e\+|e-)\d+|(\d+\.\d+|\d+\.|\.\d?|\d+))?/)?.[0]
  const unitString = input.substring(valueString.length)?.trim() || ''

  // Check value if required
  if (valueString === '') {
    valueString = undefined
    value = undefined
    if (requireValue) {
       error = 'Enter a valid numerical value'
    }
  } else {
    value = Number(valueString)
  }

  // Check unit if required
  if (requireUnit) {
    if (unitString === '') {
      return {valueString, value, error: 'Unit is required'}
    }
  }

  // Try to parse with MathJS: it can extract the unit even when it is mixed
  // with numbers
  input = normalizeExpression(input)
  let unitMathJS
  try {
    unitMathJS = UnitMathJS.parse(input, {allowNoUnits: true})
  } catch (e) {
    return {valueString, error: e.message}
  }

  let unit
  unitMathJS.value = null
  try {
    unit = new Unit(unitMathJS)
  } catch (e) {
    error = e.msg
  }
  if (error) {
    return {valueString, value, unit, error}
  }

  // If unit is not required and it is dimensionless, return without new unit
  if (!requireUnit && unit.dimension() === 'dimensionless') {
    return {valueString, value}
  }

  // TODO: This check is not enough: the input may be compatible after the base
  // units are compared.
  if (dimension !== null) {
    if (!(unit.dimension(true) === dimension || unit.dimension(false) === dimension)) {
      error = `Unit "${unit.label(false)}" is incompatible with dimension "${dimension}"`
    }
  }

  return {value, valueString, unit, error}
}
