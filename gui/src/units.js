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
import { parse, SymbolNode } from 'mathjs'
import { atom, useRecoilValue } from 'recoil'
import { scale, add } from './utils'
import { conversionMap, unitMap, unitSystems } from './unitsData'
import { isNil, isString } from 'lodash'

// Setup a unit system: by default use SI units, unless explicitly overridden
// with something else.
let defaults = {}
for (const dimension in conversionMap) {
  const info = conversionMap[dimension]
  defaults[dimension] = info.units[0]
}
const override = {
  'length': 'angstrom',
  'energy': 'electron_volt',
  'pressure': 'gigapascal',
  'angle': 'degree',
  'system': 'custom'
}
defaults = {...defaults, ...override}
export const unitsState = atom({
  key: 'units',
  default: defaults
})

const abbreviationMap = {}
for (const [key, value] of Object.entries(unitMap)) {
  abbreviationMap[value.abbreviation] = key
}

/**
 * Convenience hook for using the currently set units.
 * @returns Object containing the currently set units for each dimension (e.g.
 * {energy: 'joule'})
 */
export const useUnits = () => {
  return useRecoilValue(unitsState)
}

/**
 * Convenience function for getting the dimension of the given unit.
 *
 * @param {string||Unit} unit The unit definition for which the dimension is retrieved.
 * @returns The metainfo dimension as a string, undefined if no dimension is
 * specified.
 */
export function getDimension(unit) {
  if (unit instanceof Unit) {
    unit = unit.unitDef.name
  }
  return unitMap[unit]?.dimension
}

/**
 * Get the information related to the given unit.
 * @param {string} name Name of the unit: either the full name or an
 * abbreviation.
 *
 * @returns Information about this unit.
 */
export function getUnitByName(name) {
  // First try the full unit name
  let unitInfo = unitMap[name]
  // Try the abbreviation
  if (unitInfo === undefined) {
    const fullName = abbreviationMap[name]
    unitInfo = fullName && unitMap[fullName]
  }
  if (isNil(unitInfo)) {
    throw Error(`Unknown unit: ${name}`)
  }
  return unitInfo
}

/**
 * The the unit corresponding to the given unit in the given system.
 * @param {string} name Name of the dimension
 * abbreviation.
 *
 * @returns Information about the matching unit.
 */
export function getUnitByDimension(dimension, system) {
  if (dimension === 'dimensionless') {
    return unitMap['dimensionless']
  }

  const unitName = system[dimension]
  if (isNil(unitName)) {
    throw Error(`Could not a find unit with dimension ${dimension}" using the given unit system.`)
  }
  return getUnitByName(unitName)
}

/**
 * Helper class for persisting unit information.
 */
export class Unit {
  constructor(unit) {
    if (isString(unit)) {
      this.unitDef = this.getUnitExpression(unit)
    } else {
      this.unitDef = unit.cloneDeep()
    }
    this.label = this.getLabel(this.unitDef)
  }

  /**
   * Function for converting the value of this Unit to the SI unit system.
   *
   * @returns A new Unit instance in the SI unit system.
   */
  toSI() {
    return this.toSystem(unitSystems['SI'].units)
  }

  /**
   * Returns this unit in the given unit system.
   * @param {object} system Unit system for performing conversions.
   * @returns A new Unit instance.
   */
  toSystem(system) {
    // Create the math.js node tree from the unit definition. Unique
    // abbreviations are translated to full unit names. If a unit system is
    // given, the units are also converted here.
    let newNode = this.unitDef.cloneDeep()
    newNode = newNode.transform((node, path, parent) => {
      if (node.isSymbolNode) {
        // Get the target unit if conversion is required.
        const unitFromInfo = getUnitByName(node.name)
        const unitTo = getUnitByDimension(unitFromInfo.dimension, system)
        node.name = unitTo.name
      }
      return node
    })
    const label = this.getLabel(newNode)
    return new Unit(newNode, label)
  }

  /**
   *
   */
  getLabel(rootNode) {
    const copy = rootNode.cloneDeep()
    const newRoot = copy.transform((node, path, parent) => {
      if (node.isSymbolNode) {
        const unitFromInfo = getUnitByName(node.name)
        if (unitFromInfo !== undefined) {
          let unitToInfo = unitFromInfo
          node.name = unitToInfo.abbreviation
        }
      }
      return node
    })
    return newRoot.toString()
  }

  /**
   * Constructs a math.js expression tree for the unit in the given system. If
   * the conversion is invalid, raises an exception.
  */
  getUnitExpression(unit) {
    // Modify unit definition to comply with math.js evaluation
    let rootNode
    const from = unit.replace('**', '^')
    rootNode = parse(from)

    // Create the math.js node tree from the unit definition. Unique
    // abbreviations are translated to full unit names. If a unit system is
    // given, the units are also converted here.
    rootNode = rootNode.transform((node, path, parent) => {
      if (node.isSymbolNode) {
        // Get the unit information
        const unitFromInfo = getUnitByName(node.name)
        node.name = unitFromInfo?.name || node.name
      }
      return node
    })

    return rootNode
  }

  /**
   * Returns the system in which this unit is represented in.
   * @returns Object representing the unit system.
   */
  getUnitSystem() {
    const system = {}
    this.unitDef.transform((node, path, parent) => {
      if (node.isSymbolNode) {
        const unitFromInfo = getUnitByName(node.name)
        const dimension = unitFromInfo.dimension
        system[dimension] = unitFromInfo?.name || node.name
      }
      return node
    })
    return system
  }
}

/**
 * Helper class for persisting unit information associated with a numeric value.
 */
export class Quantity {
  constructor(value, unit) {
    this.value = value
    this.unit = unit instanceof Unit ? unit : new Unit(unit)
  }
  /**
   * Function for converting the value of this quantity to another unit.
   *
   * @param {str | Unit} unit The target unit
   * @returns A new Quantity expressed in the given units.
   */
  to(unit) {
    const newUnit = unit instanceof Unit ? unit : new Unit(unit)
    const system = newUnit.getUnitSystem()
    return this.toSystem(system)
  }
  /**
   * Function for converting the value of this Quantity to the SI unit system.
   *
   * @returns A new Quantity instance in the SI unit system.
   */
  toSI() {
    return this.toSystem(unitSystems['SI'].units)
  }
  /**
   * Function for converting the value of this quantity to another unit system.
   *
   * @param {object} system The target unit system
   * @returns A new Unit instance in the given system.
   */
  toSystem(system) {
    let value = this.value
    let unit = this.unit
    const unitDef = unit.unitDef

    // Function for getting the converted values
    function convert(value) {
      // Conversions for non-multiplicative offset units requires a special
      // handling. Here we take an automated approach: simple conversions between
      // two units takes the offset into account normally, but when combined with
      // other units (e.g.  with multiplication, division, power, etc.) the offset
      // is ignored and the 'delta' version of the units are used (see
      // https://pint.readthedocs.io/en/0.10.1/nonmult.html)
      let constant
      if (unitDef instanceof SymbolNode) {
        if (unitDef.name === 'dimensionless') {
          return value
        }
        const unitFrom = unitDef.name
        const dimension = unitMap[unitFrom].dimension
        const unitTo = system[dimension]
        constant = conversionMap[dimension]?.constants?.[unitFrom]?.[unitTo]
      }

      // Gather all units present
      const variables = new Set()
      unitDef.traverse((node, path, parent) => {
        if (node.isSymbolNode) {
          variables.add(node.name)
        }
      })

      // Check if conversion is required.
      let isConverted = true
      for (const unit of variables) {
        const dimension = unitMap[unit].dimension
        const unitTo = system[dimension]
        isConverted = unit === unitTo
        if (!isConverted) {
          break
        }
      }
      if (isConverted) {
        return value
      }

      // Gather conversion values for each present unit
      const scope = {}
      for (const unitFrom of variables) {
        const dimension = unitMap[unitFrom].dimension
        const unitTo = system[dimension]
        const multiplier = unitFrom === unitTo
          ? 1
          : conversionMap[dimension].multipliers[unitFrom][unitTo]
        scope[unitFrom] = multiplier
      }

      // Compute the scaling factor by evaluating the unit definition with the
      // SI units converted to target system
      const code = unitDef.compile()
      const factor = code.evaluate(scope)

      // Scale values to new units
      let newValues = scale(value, factor)

      // Offset if needed
      if (constant) {
        newValues = add(newValues, constant)
      }
      return newValues
    }

    return new Quantity(
      convert(value),
      new Unit(unitDef).toSystem(system)
    )
  }
}

/**
 * Convenience function for converting values to another unit system. Primarily
 * you can work with just the Unit and Quantity -classes, but sometimes it is
 * cleaner to call this function directly instead.
 *
 * @param {*} value The values to convert. Can be scalar or multi-dimensional.
 * @param {*} unit Current unit as a string or an Unit object
 * @param {object} targetUnit The target unit
 * @returns
 */
export function convertUnit(value, unit, targetUnit) {
  if (value === undefined || isNaN(value) || targetUnit === undefined || unit === undefined) return undefined
  const quantity = new Quantity(value, unit)
  const converted = quantity.to(targetUnit)
  return converted.value
}
