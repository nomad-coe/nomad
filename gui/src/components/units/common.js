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

import {memoize, has, isNil} from 'lodash'
import {Unit} from './Unit'
import { Unit as UnitMathJS } from 'mathjs'

export const deltaPrefixes = ['delta_', 'Δ']

/**
 * Convenience function for parsing value and unit information from a string.
 * Can parse values, units or both at the same time.
 *
 * @param {string} input The input string to parse
 * @param {object} options Parsing options. These include:
 *   - dimension: Dimension for the unit. Note that you should use the base
 *       dimensions which you can get e.g. with .dimension(true).
 *   - requireValue: Whether an explicit numeric value is required at the start
  *      of the input.
 *   - requireUnit: Whether an explicit unit in the input is required.
 * @returns Object containing the following properties, if available:
 *  - value: Numerical value as a number
 *  - valueString: The original number input as a string. Note that this can only return
 *      the number when it is used as a prefix, and does not work with numbers
 *      that are part of a complex expression, e.g. 300 eV / 1000 K.
 *  - unit: Unit instance
 *  - error: Error messsage
 */
export function parse(input, options) {
  options = {
    dimension: null,
    requireValue: false,
    requireUnit: false,
    ...options
  } || {}

  const result = {}

  try {
    const parseResults = parseInternal(input, options)
    result.value = isNil(parseResults.value) ? undefined : parseResults.value
    result.valueString = parseResults.valueString || undefined
    if (parseResults.unit?.units?.length) {
      result.unit = new Unit(parseResults.unit)
    }
  } catch (e) {
    result.error = e.message
    return result
  }

  // TODO: This check is not enough: the input may be compatible after the base
  // units are compared.
  if (options.dimension !== null && result.unit) {
    if (!(result.unit.dimension(true) === options.dimension || result.unit.dimension(false) === options.dimension)) {
      result.error = `Unit "${result.unit.label(false)}" is incompatible with dimension "${options.dimension}".`
    }
  }

  return result
}

/**
 * Parse a string into an optional value and a MathJS Unit definition.
 *
 * Throws an exception if the provided string does not contain a valid unit or
 * cannot be parsed.
 *
 * @memberof Unit
 * @param {string} str A string like "5.2 inch", "4e2 cm/s^2"
 * @param {object} options Parsing options. These include:
 *   - requireValue: Whether an explicit numeric value is required at the start
  *      of the input.
 *   - requireUnit: Whether an explicit unit in the input is required.
 * @returns Object containing the following properties, if available:
 *  - value: Numerical value as a number
 *  - valueString: The original number input as a string. Note that this can only return
 *      the number when it is used as a prefix, and does not work with numbers
 *      that are part of a complex expression, e.g. 300 eV / 1000 K.
 *  - unit: Unit instance
 */
export function parseInternal(str, options) {
  // Replace ** with ^
  str = str.replace(/\*\*/g, '^')

  options = options || {}
  const text = str
  let index = -1
  let c = ''

  function skipWhitespace() {
    while (c === ' ' || c === '\t') {
      next()
    }
  }

  function isDigitDot(c) {
    return ((c >= '0' && c <= '9') || c === '.')
  }

  function isDigit(c) {
    return ((c >= '0' && c <= '9'))
  }

  function next() {
    index++
    c = text.charAt(index)
  }

  function revert(oldIndex) {
    index = oldIndex
    c = text.charAt(index)
  }

  function parseUnit() {
    let unitName = ''

    // Alphanumeric characters only; matches [a-zA-Z0-9]
    while (isDigit(c) || UnitMathJS.isValidAlpha(c)) {
      unitName += c
      next()
    }

    // Must begin with [a-zA-Z]
    const firstC = unitName.charAt(0)
    if (UnitMathJS.isValidAlpha(firstC)) {
      return unitName
    } else {
      return null
    }
  }

  function parseCharacter(toFind) {
    if (c === toFind) {
      next()
      return toFind
    } else {
      return null
    }
  }

  function parseNumber() {
    let number = ''
    const oldIndex = index

    if (c === '+') {
      next()
    } else if (c === '-') {
      number += c
      next()
    }

    if (!isDigitDot(c)) {
      // a + or - must be followed by a digit
      revert(oldIndex)
      return null
    }

    // get number, can have a single dot
    if (c === '.') {
      number += c
      next()
      if (!isDigit(c)) {
        // this is no legal number, it is just a dot
        revert(oldIndex)
        return null
      }
    } else {
      while (isDigit(c)) {
        number += c
        next()
      }
      if (c === '.') {
        number += c
        next()
      }
    }
    while (isDigit(c)) {
      number += c
      next()
    }

    // check for exponential notation like "2.3e-4" or "1.23e50"
    if (c === 'E' || c === 'e') {
      // The grammar branches here. This could either be part of an exponent or the start of a unit that begins with the letter e, such as "4exabytes"

      let tentativeNumber = ''
      const tentativeIndex = index

      tentativeNumber += c
      next()

      if (c === '+' || c === '-') {
        tentativeNumber += c
        next()
      }

      // Scientific notation MUST be followed by an exponent (otherwise we assume it is not scientific notation)
      if (!isDigit(c)) {
        // The e or E must belong to something else, so return the number without the e or E.
        revert(tentativeIndex)
        return number
      }

      // We can now safely say that this is scientific notation.
      number = number + tentativeNumber
      while (isDigit(c)) {
        number += c
        next()
      }
    }

    return number
  }

  if (typeof text !== 'string') {
    throw new TypeError('Invalid argument in Unit.parse, string expected')
  }

  const unit = new UnitMathJS()
  unit.units = []

  let powerMultiplierCurrent = 1
  let expectingUnit = false

  // A unit should follow this pattern:
  // [number] ...[ [*/] unit[^number] ]
  // unit[^number] ... [ [*/] unit[^number] ]

  // Rules:
  // number is any floating point number.
  // unit is any alphanumeric string beginning with an alpha. Units with names like e3 should be avoided because they look like the exponent of a floating point number!
  // The string may optionally begin with a number.
  // Each unit may optionally be followed by ^number.
  // Whitespace or a forward slash is recommended between consecutive units, although the following technically is parseable:
  //   2m^2kg/s^2
  // it is not good form. If a unit starts with e, then it could be confused as a floating point number:
  //   4erg

  next()
  skipWhitespace()

  // Optional number at the start of the string
  const valueString = parseNumber()
  let value = null
  if (valueString) {
    value = parseFloat(valueString)
    skipWhitespace() // Whitespace is not required here

    // handle multiplication or division right after the value, like '1/s'
    if (parseCharacter('*')) {
      powerMultiplierCurrent = 1
      expectingUnit = true
    } else if (parseCharacter('/')) {
      powerMultiplierCurrent = -1
      expectingUnit = true
    }
  } else if (options.requireValue) {
    throw new SyntaxError('Enter a valid numerical value')
  }

  // Stack to keep track of powerMultipliers applied to each parentheses group
  const powerMultiplierStack = []

  // Running product of all elements in powerMultiplierStack
  let powerMultiplierStackProduct = 1

  while (true) {
    skipWhitespace()

    // Check for and consume opening parentheses, pushing powerMultiplierCurrent to the stack
    // A '(' will always appear directly before a unit.
    while (c === '(') {
      powerMultiplierStack.push(powerMultiplierCurrent)
      powerMultiplierStackProduct *= powerMultiplierCurrent
      powerMultiplierCurrent = 1
      next()
      skipWhitespace()
    }

    // Is there something here?
    let uStr
    if (c) {
      const oldC = c
      uStr = parseUnit()
      if (uStr === null) {
        throw new SyntaxError('Unexpected "' + oldC + '" in "' + text + '" at index ' + index.toString())
      }
    } else {
      // End of input.
      break
    }

    // Verify the unit exists and get the prefix (if any)
    const res = findUnit(uStr)
    if (res === null) {
      // Unit not found.
      throw new SyntaxError('Unit "' + uStr + '" not found.')
    }

    let power = powerMultiplierCurrent * powerMultiplierStackProduct
    // Is there a "^ number"?
    skipWhitespace()
    if (parseCharacter('^')) {
      skipWhitespace()
      const p = parseNumber()
      if (p === null) {
        // No valid number found for the power!
        throw new SyntaxError('In "' + str + '", "^" must be followed by a floating-point number')
      }
      power *= p
    }

    // Add the unit to the list
    unit.units.push({
      unit: res.unit,
      prefix: res.prefix,
      delta: res.delta,
      power
    })
    for (let i = 0; i < UnitMathJS.BASE_DIMENSIONS.length; i++) {
      unit.dimensions[i] += (res.unit.dimensions[i] || 0) * power
    }

    // Check for and consume closing parentheses, popping from the stack.
    // A ')' will always follow a unit.
    skipWhitespace()
    while (c === ')') {
      if (powerMultiplierStack.length === 0) {
        throw new SyntaxError('Unmatched ")" in "' + text + '" at index ' + index.toString())
      }
      powerMultiplierStackProduct /= powerMultiplierStack.pop()
      next()
      skipWhitespace()
    }

    // "*" and "/" should mean we are expecting something to come next.
    // Is there a forward slash? If so, negate powerMultiplierCurrent. The next unit or paren group is in the denominator.
    expectingUnit = false

    if (parseCharacter('*')) {
      // explicit multiplication
      powerMultiplierCurrent = 1
      expectingUnit = true
    } else if (parseCharacter('/')) {
      // division
      powerMultiplierCurrent = -1
      expectingUnit = true
    } else {
      // implicit multiplication
      powerMultiplierCurrent = 1
    }

    // Replace the unit into the auto unit system
    if (res.unit.base) {
      const baseDim = res.unit.base.key
      UnitMathJS.UNIT_SYSTEMS.auto[baseDim] = {
        unit: res.unit,
        prefix: res.prefix
      }
    }
  }

  // Has the string been entirely consumed?
  skipWhitespace()
  if (c) {
    throw new SyntaxError('Could not parse: "' + str + '"')
  }

  // Is there a trailing slash?
  if (expectingUnit) {
    throw new SyntaxError('Trailing characters: "' + str + '"')
  }

  // Is the parentheses stack empty?
  if (powerMultiplierStack.length !== 0) {
    throw new SyntaxError('Unmatched "(" in "' + text + '"')
  }

  // Are there any units at all?
  if (unit.units.length === 0 && options.requireUnit) {
    throw new SyntaxError('Unit is required')
  }

  return {value, valueString, unit}
}

/**
 * Find a unit from a string
 *
 * @param {string} str A string like 'cm' or 'inch'
 * @returns {Object | null} When found, an object with fields unit and
*    prefix is returned. Else, null is returned.
 */
const findUnit = memoize((str) => {
  // First, match units names exactly. For example, a user could define 'mm' as
  // 10^-4 m, which is silly, but then we would want 'mm' to match the
  // user-defined unit.
  if (has(UnitMathJS.UNITS, str)) {
    const unit = UnitMathJS.UNITS[str]
    const prefix = unit.prefixes['']
    return { unit, prefix, delta: false }
  }

  for (const name in UnitMathJS.UNITS) {
    if (has(UnitMathJS.UNITS, name)) {
      if (str.endsWith(name)) {
        const unit = UnitMathJS.UNITS[name]
        const prefixLen = (str.length - name.length)
        let prefixName = str.substring(0, prefixLen)
        let delta = false
        for (const deltaPrefix of deltaPrefixes) {
          if (prefixName.startsWith(deltaPrefix)) {
            prefixName = prefixName.substring(deltaPrefix.length)
            delta = true
            break
          }
        }
        const prefix = has(unit.prefixes, prefixName)
          ? unit.prefixes[prefixName]
          : undefined
        if (prefix !== undefined) {
          return { unit, prefix, delta }
        }
      }
    }
  }

  return null
})
