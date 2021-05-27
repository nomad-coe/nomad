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

// Set up a unit system: by default use SI units, unless explicitly overridden
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
 * @returns Object containing the currently set unit system.
 */
export const useUnits = () => {
  return useRecoilValue(unitsState)
}

/**
 * Helper class for persisting unit information.
 */
export class Unit {
  constructor(unit, system) {
    this.unitDef = getUnitDefinition(unit, system)
  }
  label() {
    const copy = this.unitDef.cloneDeep()
    const newRoot = copy.transform((node, path, parent) => {
      if (node.isSymbolNode) {
        node.name = unitMap[node.name].abbreviation
      }
      return node
    })
    return newRoot.toString()
  }
}

/**
 * Helper class for persisting unit information for a numeric value.
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
 *
 * @param {*} value
 * @param {*} unitDef
 * @param {*} system
 * @returns
 */
export function toUnitSystem(value, unit, system) {
  // If unit is not already given as a math.js nodetree, convert it to such.
  let unitDef
  if (!(unit instanceof Unit)) {
    unitDef = getUnitDefinition(unit, system)
  } else {
    unitDef = unit.unitDef
  }
  // Temperatures require special handling due to the fact that Celsius and
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

export function convertSILabel(unit, system) {
  // Modify unit definition to comply with math.js evaluation
  const from = unit.replace('**', '^')

  // Form new unit definition string by replacing the SI units with the target
  const rootNode = parse(from)
  const newRoot = rootNode.transform((node, path, parent) => {
    if (node.isSymbolNode) {
      const unitFromInfo = unitMap[node.name]
      if (unitFromInfo !== undefined) {
        const dimension = unitFromInfo.dimension
        const unitTo = system[dimension]
        const unitToInfo = unitMap[unitTo]
        const label = unitToInfo.abbreviation
        node.name = label
      }
    }
    return node
  })
  return newRoot.toString()
}

export function getUnitDefinition(unit, system) {
  // Modify unit definition to comply with math.js evaluation
  const from = unit.replace('**', '^')

  // Form new unit definition string by replacing the SI units with the target
  const rootNode = parse(from)
  const newRoot = rootNode.transform((node, path, parent) => {
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
  return newRoot
}

/**
 * Used to convert a unit definition into an human-friendly abbreviated form.
 *
 * @param {string} unit The unit definition to convert. Use the standardized
 * unit names used in unit.js
 * @returns A human-friendly abbreviated form of the unit.
 */
export function getUnitLabel(unitDefinition) {
  const copy = unitDefinition.cloneDeep()
  const newRoot = copy.transform((node, path, parent) => {
    if (node.isSymbolNode) {
      node.name = unitMap[node.name].abbreviation
    }
    return node
  })
  return newRoot.toString()
}

/**
 * Used to convert numeric values from SI units to the given unit system. Works
 * on n-dimensional arrays.
 *
 * @param {*} value The values to convert. Can be a scalar or an n-dimensional
 * array.
 * @param {string} unit Original SI unit definition. Can be any algebraic
 * combination of SI units, e.g. "1 / meter^2". The unit names should follow
 * the definitions provided in the file units.js that is generated by the NOMAD
 * CLI.
 * @param {*} system Target unit system. A Javascript object where each
 * physical quantity (e.g. "length") acts as a key that corresponds to a target
 * unit (e.g. "angstrom"). The unit names should follow the definitions
 * provided in the file units.js that is generated by the NOMAD CLI.
 * array.
 * @param {string} units Whether to return the new unit as a string together
 * with the new values.
 *
 * @return {*} A copy of the original data with units converted.
 */
export function convertSI(value, unit, system, units = true) {
  // Modify unit definition to comply with math.js evaluation
  const from = unit.replace('**', '^')

  // Temperatures require special handling due to the fact that Celsius and
  // Fahrenheit are not absolute units and are non-multiplicative. Two kinds of
  // temperature conversions are supported: ones with a single temperature unit
  // and ones where temperature is used as a part of an expression. If a single
  // temperature unit is specified, they are converted normally taking the
  // offset into account. If they are used as a part of an expression, they are
  // interpreted as ranges and the offset is ignored.
  if (from === 'kelvin') {
    const unitTo = system.units['temperature']
    const multiplier = conversionMap['temperature'].multipliers['kelvin'][unitTo]
    const constant = conversionMap['temperature'].constants['kelvin'][unitTo]
    const label = unitMap['kelvin'].label
    let newValues = value
    if (multiplier !== 1) {
      newValues = scale(newValues, multiplier)
    }
    if (constant !== undefined) {
      newValues = add(newValues, constant)
    }
    if (units) {
      return [newValues, label]
    }
    return newValues
  }

  // Gather all units present
  const variables = new Set()
  const rootNode = parse(from)
  rootNode.traverse((node, path, parent) => {
    if (node.isSymbolNode) {
      variables.add(node.name)
    }
  })

  // Check if conversion is required. The unit definition string is standardized
  // even if no conversion took place.
  let isSI = true
  for (const unit of variables) {
    const dimension = unitMap[unit].dimension
    const unitSI = unitSystems['SI'].units[dimension]
    isSI = unit === unitSI
    if (!isSI) {
      break
    }
  }
  if (isSI) {
    if (units) {
      const newUnit = convertSILabel(from, system)
      return [value, newUnit]
    }
    return value
  }

  // Gather conversion values for each present SI unit
  const scope = {}
  for (const unitFrom of variables) {
    const dimension = unitMap[unitFrom].dimension
    const unitTo = system.units[dimension]
    scope[unitFrom] = conversionMap[dimension].multipliers[unitFrom][unitTo]
  }

  // Compute the scaling factor by evaluating the unit definition with the
  // SI units converted to target system
  const code = rootNode.compile()
  const factor = code.evaluate(scope)

  // Scale values to new units
  let newValues = scale(value, factor)

  // Form new unit definition string by replacing the SI units with the target
  if (units) {
    const newUnit = convertSILabel(from, system)
    return [newValues, newUnit]
  }
  return newValues
}
