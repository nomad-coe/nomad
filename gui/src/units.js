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
  'system': 'custom'
}
defaults = {...defaults, ...override}
export const unitsState = atom({
  key: 'units',
  default: defaults
})

/**
 * Convenience hook for using the currently set units.
 * @returns Object containing the currently set units for each dimension (e.g.
 * {energy: 'joule'})
 */
export const useUnits = () => {
  return useRecoilValue(unitsState)
}

/**
 * Helper class for persisting unit information.
 */
export class Unit {
  constructor(unit, system = undefined) {
    this.unitDef = this.getUnitDefinition(unit, system)
  }
  getUnitDefinition(unit, system = undefined) {
    // Modify unit definition to comply with math.js evaluation
    const from = unit.replace('**', '^')

    // Create the math.js node tree from the unit definition. If a unit system
    // is given, the units are also converted here.
    let rootNode = parse(from)
    if (system) {
      rootNode = rootNode.transform((node, path, parent) => {
        if (node.isSymbolNode) {
          const unitFromInfo = unitMap[node.name]
          if (unitFromInfo !== undefined) {
            const dimension = unitFromInfo.dimension
            const unitTo = system[dimension]
            node.name = unitTo
          }
        }
        return node
      })
    }
    return rootNode
  }
  /**
   * Returns the unit definition as a human-readable string. If a unit system is
   * provided, performs the conversions to that system on-the-fly.
   * @param {object} system Optional unit system for performing conversions.
   * @returns The unit definition as a string.
   */
  label(system = undefined) {
    const copy = this.unitDef.cloneDeep()
    const newRoot = copy.transform((node, path, parent) => {
      if (node.isSymbolNode) {
        const unitFromInfo = unitMap[node.name]
        if (unitFromInfo !== undefined) {
          let def
          if (system) {
            const dimension = unitFromInfo.dimension
            const unitTo = system[dimension]
            def = unitMap[unitTo]
          } else {
            def = unitFromInfo
          }
          node.name = def.abbreviation
        }
      }
      return node
    })
    return newRoot.toString()
  }
}

/**
 * Helper class for persisting unit information associated with a numeric value.
 */
export class Quantity {
  constructor(value, unit) {
    this.value = value
    this.unit = unit
  }
  toSI() {
    return toSI(this.value, this.unit)
  }
  toSystem(system) {
    return toUnitSystem(this.value, this.unit, system)
  }
}

/**
 * Convenience function for converting values to another unit system. Primarily
 * you can work with just the Unit and Quantity -classes, but sometimes it is
 * cleaner to call this function directly instead.
 *
 * @param {*} value
 * @param {*} unitDef
 * @param {*} system
 * @returns
 */
export function toUnitSystem(value, unit, system) {
  // If value given as Quantity, extract unit and value from it
  if (value instanceof Quantity) {
    value = value.value
    unit = value.unit
  }

  // If unit is given as a string, create a Unit from for it.
  if (!(unit instanceof Unit)) {
    unit = new Unit(unit)
  }
  const unitDef = unit.unitDef
  // Temperatures require special handling due to the fact that e.g. Celsius and
  // Fahrenheit are not absolute units and are non-multiplicative. Two kinds of
  // temperature conversions are supported: ones with a single temperature unit
  // and ones where temperature is used as a part of an expression. If a single
  // temperature unit is specified, they are converted normally taking the
  // offset into account. If they are used as a part of an expression, they are
  // interpreted as ranges and the offset is ignored.
  if (unitDef instanceof SymbolNode) {
    const unitFrom = unitDef.name
    if (unitMap[unitFrom].dimension === 'temperature') {
      const unitTo = system['temperature']
      const multiplier = conversionMap['temperature'].multipliers[unitFrom][unitTo]
      const constant = conversionMap['temperature'].constants[unitFrom][unitTo]
      let newValues = value
      if (multiplier !== 1) {
        newValues = scale(newValues, multiplier)
      }
      if (constant !== undefined) {
        newValues = add(newValues, constant)
      }
      return newValues
    }
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
    scope[unitFrom] = conversionMap[dimension].multipliers[unitFrom][unitTo]
  }

  // Compute the scaling factor by evaluating the unit definition with the
  // SI units converted to target system
  const code = unitDef.compile()
  const factor = code.evaluate(scope)

  // Scale values to new units
  let newValues = scale(value, factor)
  return newValues
}

/**
 *
 * @param {*} value
 * @param {*} unitDef
 * @param {*} system
 * @returns
 */
export function toSI(value, unit) {
  return toUnitSystem(value, unit, unitSystems['SI'].units)
}
