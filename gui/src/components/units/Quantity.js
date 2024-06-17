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
      throw Error('Please provide the value for a Quantity as a number, or as a multidimensional array of numbers.')
    }

    // This attribute stores the quantity value in 'normalized' form that is
    // given in the base units (=SI). This value should only be determined once
    // during initialization and all calls to value() will then lazily determine
    // the value in the currently set units. This avoids 'drift' in the value
    // caused by several consecutive changes of the units.
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
  normalize(values) {
    // Pre-calculate coefficients for currently set units. This speeds up
    // conversion for large arrays.
    const unit = this.unit.mathjsUnit
    const ignoreOffset = unit._isDerived()
    const coefficients = this.conversion_coefficients()

    return mapDeep(values, (value) => {
      if (value === null || value === undefined) {
        return value
      }

      let result = value
      for (let i = 0; i < unit.units.length; i++) {
        const unitDef = unit.units[i]
        const unitOffset = unitDef.unit.offset
        const variable = (ignoreOffset || unitDef.delta)
          ? result
          : result + unitOffset
        result = variable * coefficients[i]
      }
      return result
    })
  }

  /**
   * Convert value from base units to currently set units.
   * @param {n-dimensional array} value Value in base units.
   * @returns Value in currently set units.
   */
  denormalize(values) {
    // Pre-calculate coefficients for currently set units. This speeds up
    // conversion for large arrays.
    const unit = this.unit.mathjsUnit
    const ignoreOffset = unit._isDerived()
    const coefficients = this.conversion_coefficients()

    return mapDeep(values, (value) => {
      if (value === null || value === undefined) {
        return value
      }

      let result = value
      for (let i = 0; i < unit.units.length; i++) {
        const unitDef = unit.units[i]
        const unitOffset = unitDef.unit.offset
        result = (ignoreOffset || unitDef.delta)
          ? result / coefficients[i]
          : result / coefficients[i] - unitOffset
      }
      return result
    })
  }

  /**
   * Returns a set of conversion coefficients based on the currently set units.
   * @returns Array of conversion coefficients, one for each unit that is present.
   */
  conversion_coefficients() {
    const unit = this.unit.mathjsUnit
    const coefficients = []
    for (let i = 0; i < unit.units.length; i++) {
      const unitDef = unit.units[i]
      const unitValue = unitDef.unit.value
      const unitPrefixValue = unitDef.prefix.value
      const unitPower = unitDef.power
      coefficients.push((Math.pow(unitValue * unitPrefixValue, unitPower)))
    }
    return coefficients
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
