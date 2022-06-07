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
import { isNil, startCase, toLower, has, cloneDeep, isString, isNumber, isArray } from 'lodash'
import { deepMap } from './utils'
import { Unit as UnitMathJS, createUnit } from 'mathjs'
import { atom, useRecoilValue } from 'recoil'
import { unitList, prefixes } from './unitsData'

// Delete all units and prefixes that come by default with Math.js. This way
// they cannot be intermixed with the NOMAD units. Notice that we have to clear
// them in place: they are defined as const.
Object.keys(UnitMathJS.UNITS).forEach(name => UnitMathJS.deleteUnit(name))
UnitMathJS.BASE_DIMENSIONS.splice(0, UnitMathJS.BASE_DIMENSIONS.length)
Object.getOwnPropertyNames(UnitMathJS.BASE_UNITS).forEach(function(prop) {
  delete UnitMathJS.BASE_UNITS[prop]
})
Object.getOwnPropertyNames(UnitMathJS.PREFIXES).forEach(function(prop) {
  delete UnitMathJS.PREFIXES[prop]
})
UnitMathJS.PREFIXES.NONE = {'': { name: '', value: 1, scientific: true }}
UnitMathJS.PREFIXES.PINT = prefixes

// Customize the unit parsing to allow certain special symbols
const isAlphaOriginal = UnitMathJS.isValidAlpha
const isSpecialChar = function(c) {
  const specialChars = new Set(['_', 'Å', 'Å', 'å', '°', 'µ', 'ö', 'é', '∞'])
  return specialChars.has(c)
}
const isGreekLetter = function(c) {
  const charCode = c.charCodeAt(0)
  return (charCode > 912 && charCode < 970)
}
UnitMathJS.isValidAlpha = function(c) {
  return isAlphaOriginal(c) || isSpecialChar(c) || isGreekLetter(c)
}

// Create MathJS unit definitions from the data exported by 'nomad dev units'
const unitToAbbreviationMap = {}
const unitDefinitions = {}
for (let def of unitList) {
  const name = def.name
  def = {
    ...def,
    baseName: def.dimension,
    prefixes: 'pint'
  }
  unitDefinitions[name] = def
  // Register abbreviations
  if (def.abbreviation) {
    unitToAbbreviationMap[name] = def.abbreviation
    if (def.aliases) {
      for (const alias of def.aliases) {
        unitToAbbreviationMap[alias] = def.abbreviation
      }
    }
  }
}
createUnit(unitDefinitions, {override: true})

export const unitMap = Object.fromEntries(unitList.map(x => [x.name, x]))

// Define the unit systems
const SIUnits = {
  // Base units
  dimensionless: {name: 'dimensionless', fixed: false},
  length: {name: 'meter', fixed: false},
  mass: {name: 'kilogram', fixed: false},
  time: {name: 'second', fixed: false},
  current: {name: 'ampere', fixed: false},
  temperature: {name: 'kelvin', fixed: false},
  luminosity: {name: 'candela', fixed: false},
  luminous_flux: {name: 'lumen', fixed: false},
  substance: {name: 'mole', fixed: false},
  angle: {name: 'radian', fixed: false},
  information: {name: 'bit', fixed: false},
  // Derived units with specific name
  force: {name: 'newton', fixed: false},
  energy: {name: 'joule', fixed: false},
  power: {name: 'watt', fixed: false},
  pressure: {name: 'pascal', fixed: false},
  charge: {name: 'coulomb', fixed: false},
  resistance: {name: 'ohm', fixed: false},
  conductance: {name: 'siemens', fixed: false},
  inductance: {name: 'henry', fixed: false},
  magnetic_flux: {name: 'weber', fixed: false},
  magnetic_field: {name: 'tesla', fixed: false},
  frequency: {name: 'hertz', fixed: false},
  luminance: {name: 'nit', fixed: false},
  illuminance: {name: 'lux', fixed: false},
  electric_potential: {name: 'volt', fixed: false},
  capacitance: {name: 'farad', fixed: false},
  activity: {name: 'katal', fixed: false}
  // TODO: Derived units without a specific name. Cannot be registered as such
  // as they don't have labels or separate definitions.
  // area: {name: 'm^2', fixed: false},
  // volume: {name: 'm^2', fixed: false},
  // wavenumber: {name: '1 / m', fixed: false},
  // speed: {name: 'm/s', fixed: false},
  // acceleration: {name: 'm / s^2', fixed: false},
  // density: {name: 'kg / m^3', fixed: false},
  // viscosity: {name: 'Pa / s', fixed: false},
  // kinematic_viscosity: {name: 'm^2 / s', fixed: false},
  // fluidity: {name: '1 / Pa / s', fixed: false},
  // concentration: {name: 'mol / m^3', fixed: false},
  // entropy: {name: 'J / K', fixed: false},
  // molar_entropy: {name: 'J / K / mole', fixed: false},
  // electric_field: {name: 'V / m', fixed: false},
  // intensity: {name: 'W / m^2', fixed: false},
  // electric_dipole: {name: 'C m', fixed: false},
  // electric_quadrupole: {name: 'C m^2', fixed: false},
  // magnetic_dipole: {name: 'A m^2', fixed: false}
}
const SIUnitsFixed = cloneDeep(SIUnits)
Object.values(SIUnitsFixed).forEach(value => { value.fixed = true })

export const unitSystems = {
  Custom: {
    label: 'Custom',
    description: 'Custom unit setup',
    units: {
      ...SIUnits,
      length: {name: 'meter', fixed: false},
      time: {name: 'second', fixed: false},
      energy: {name: 'electron_volt', fixed: false},
      pressure: {name: 'pascal', fixed: false},
      angle: {name: 'degree', fixed: false}
    }
  },
  SI: {
    label: 'SI',
    description: 'International System of Units (SI)',
    units: SIUnitsFixed
  },
  AU: {
    label: 'AU',
    description: 'Hartree atomic units',
    units: {
      ...SIUnits,
      time: {name: 'atomic_unit_of_time', fixed: true},
      length: {name: 'bohr', fixed: true},
      mass: {name: 'electron_mass', fixed: true},
      current: {name: 'atomic_unit_of_current', fixed: true},
      temperature: {name: 'atomic_unit_of_temperature', fixed: true},
      force: {name: 'atomic_unit_of_force', fixed: true},
      energy: {name: 'hartree', fixed: true},
      pressure: {name: 'atomic_unit_of_pressure', fixed: true},
      angle: {name: 'radian', fixed: true}
    }
  }
}

// Create a map of all units per dimension
export const dimensionMap = {}
for (const def of unitList) {
  const name = def.name
  const dimension = def.dimension
  if (isNil(dimension)) {
    continue
  }
  const oldInfo = dimensionMap[dimension] || {
    label: startCase(toLower(dimension.replace('_', ' ')))
  }
  const oldList = oldInfo.units || []
  oldList.push(name)
  oldInfo.units = oldList
  dimensionMap[dimension] = oldInfo
}

// Check that all units in the unit systems have been registered
for (const [systemName, system] of Object.entries(unitSystems)) {
  for (const [dimension, unit] of Object.entries(system.units)) {
    if (isNil(UnitMathJS.UNITS[unit.name])) {
      throw Error(`Unknown unit for dimension '${dimension}' found in system '${systemName}'`)
    }
  }
}

// A state containing the currently configured unit system.
export const unitsState = atom({
  key: 'units',
  default: unitSystems.Custom
})

/**
 * Convenience hook for using the currently set units.
 * @returns Object containing the currently set units for each dimension (e.g.
 * {energy: 'joule'})
 */
export const useUnits = () => {
  const unitSystem = useRecoilValue(unitsState)
  return unitSystem.units
}

/**
 * Helper class for persisting unit information.
 *
 * Builds upon the math.js Unit class system, but adds additional functionality,
 * including:
 *  - Ability to convert to any unit system given as an argument
 *  - Abbreviated labels for dense formatting
 */
export class Unit {
  /**
   * @param {str | Unit} unit Unit for the quantity.
   */
  constructor(unit) {
    if (isString(unit)) {
      unit = this.normalizeExpression(unit)
      unit = new UnitMathJS(undefined, unit)
    } else if (unit instanceof Unit) {
      unit = unit.mathjsUnit.clone()
    } else if (unit instanceof UnitMathJS) {
      unit = unit.clone()
    } else {
      throw Error('Please provide the unit as a string or as an instance of Unit.')
    }
    this.mathjsUnit = unit
    // this._labelabbreviate = undefined
    // this._label = undefined
  }

  /**
   * Normalizes the given expression into a format that can be parsed by MathJS.
   * @param {str} expression Expression
   * @returns string Expression in normalized form
   */
  normalizeExpression(expression) {
    return expression.replace('**', '^')
  }

  /**
   * Checks if the given unit has the same based dimensions as this one.
   * @param {str | Unit} unit Unit to compare to
   * @returns boolean Whether the units have the same base dimensions.
   */
  equalBase(unit) {
    if (isString(unit)) {
      unit = this.normalizeExpression(unit)
      unit = new Unit(unit)
    }
    return this.mathjsUnit.equalBase(unit.mathjsUnit)
  }

  /**
   * Used to create a human-readable description of the unit as a string.
   *
   * @param {bool} abbreviate Whether to abbreviate the label using the
   * abbreviations for each unit and prefix. If false, the original unit names
   * (as given or defined by the unit system) are used.
   * @returns A string representing the unit.
   */
  label(abbreviate = true) {
    // TODO: The label caching is disabled for now. Because Quantities are
    // stored as recoil.js atoms, they become immutable which causes problems
    // with internal state mutation.
    // if (this._labelabbreviate === abbreviate && this._label) {
    //   return this._label
    // }
    const units = this.mathjsUnit.units
    let strNum = ''
    let strDen = ''
    let nNum = 0
    let nDen = 0

    function getName(unit) {
      if (unit.base.key === 'dimensionless') {
        return ''
      }
      if (!abbreviate) {
        return unit.name
      }
      return unitToAbbreviationMap?.[unit.name] || unit.name
    }

    function getPrefix(unit, original) {
      if (!abbreviate) {
        return original
      }
      const prefixMap = {
        // SI
        deca: 'da',
        hecto: 'h',
        kilo: 'k',
        mega: 'M',
        giga: 'G',
        tera: 'T',
        peta: 'P',
        exa: 'E',
        zetta: 'Z',
        yotta: 'Y',
        deci: 'd',
        centi: 'c',
        milli: 'm',
        micro: 'u',
        nano: 'n',
        pico: 'p',
        femto: 'f',
        atto: 'a',
        zepto: 'z',
        yocto: 'y',
        // IEC
        kibi: 'Ki',
        mebi: 'Mi',
        gibi: 'Gi',
        tebi: 'Ti',
        pebi: 'Pi',
        exi: 'Ei',
        zebi: 'Zi',
        yobi: 'Yi'
      }
      return prefixMap?.[original] || original
    }

    for (let i = 0; i < units.length; i++) {
      if (units[i].power > 0) {
        nNum++
        const prefix = getPrefix(units[i].unit.name, units[i].prefix.name)
        const name = getName(units[i].unit)
        strNum += ` ${prefix}${name}`
        if (Math.abs(units[i].power - 1.0) > 1e-15) {
          strNum += '^' + units[i].power
        }
      } else if (units[i].power < 0) {
        nDen++
      }
    }

    if (nDen > 0) {
      for (let i = 0; i < units.length; i++) {
        if (units[i].power < 0) {
          const prefix = getPrefix(units[i].unit.name, units[i].prefix.name)
          const name = getName(units[i].unit)
          if (nNum > 0) {
            strDen += ` ${prefix}${name}`
            if (Math.abs(units[i].power + 1.0) > 1e-15) {
              strDen += '^' + (-units[i].power)
            }
          } else {
            strDen += ` ${prefix}${name}`
            strDen += '^' + (units[i].power)
          }
        }
      }
    }
    // Remove leading whitespace
    strNum = strNum.substr(1)
    strDen = strDen.substr(1)

    // Add parentheses for better copy/paste back into evaluate, for example, or
    // for better pretty print formatting
    if (nNum > 1 && nDen > 0) {
      strNum = '(' + strNum + ')'
    }
    if (nDen > 1 && nNum > 0) {
      strDen = '(' + strDen + ')'
    }

    let str = strNum
    if (nNum > 0 && nDen > 0) {
      str += ' / '
    }
    str += strDen

    // this._labelabbreviate = abbreviate
    // this._label = str
    return str
  }

  /**
   * Gets the dimension of this unit as a string. The order of the dimensions is
   * fixed (determined at unit registration time).
   *
   * @param {boolean} base Whether to return dimension in base units. Othwerwise
   * the original unit dimensions are used.
   * @returns The dimensionality as a string, e.g. 'time^2 energy mass^-2'
   */
  dimension(base = true) {
    const dimensions = Object.keys(UnitMathJS.BASE_UNITS)
    const dimensionMap = Object.fromEntries(dimensions.map(name => [name, 0]))

    if (base) {
      const BASE_DIMENSIONS = UnitMathJS.BASE_DIMENSIONS
      for (let i = 0; i < BASE_DIMENSIONS.length; ++i) {
        const power = this?.mathjsUnit.dimensions?.[i]
        if (power) {
          dimensionMap[BASE_DIMENSIONS[i]] += power
        }
      }
    } else {
      for (const unit of this?.mathjsUnit.units) {
        const power = unit.power
        if (power) {
          dimensionMap[unit.unit.base.key] += power
        }
      }
    }
    return Object.entries(dimensionMap)
      .filter(d => d[1] !== 0)
      .map(d => `${d[0]}${((d[1] < 0 || d[1] > 1) && `^${d[1]}`) || ''}`).join(' ')
  }

  /**
   * Function for converting to another unit.
   *
   * @param {str | Unit} unit The target unit
   * @returns A new Unit expressed in the given units.
   */
  to(unit) {
    if (isString(unit)) {
      unit = this.normalizeExpression(unit)
    } else if (unit instanceof Unit) {
      unit = unit.label()
    } else {
      throw Error('Unknown unit type. Please provide the unit as as string or as instance of Unit.')
    }
    return new Unit(this.mathjsUnit.to(unit))
  }

  /**
   * Function for converting the value of this Unit to the SI unit system.
   *
   * @returns A new Unit instance in the SI unit system.
   */
  toSI() {
    return this.toSystem(unitSystems.SI.units)
  }

  /**
   * Function for converting the value of this unit to another unit system.
   *
   * Notice that converting a unit to another unit system is not as easy as
   * conversions to a specific unit. When converting to a specific unit one can
   * simply check that the dimensions match and go ahead with the conversion.
   * With unit systems, there can be multiple alternative forms, and choosing a
   * good one is more difficult. E.g. should 'a_u_force * angstrom' be converted
   * into:
   *
   * a) N m
   * b) J
   * c) (kg m^2) / s^2
   *
   * By default this function will try to preserve the original unit dimensions
   * and not convert everything down to base units. If a derived unit is not
   * present, it will, however, attempt to convert it to the base units. Any
   * further simplication is not performed.
   *
   * @param {object} system The target unit system.
   * @returns A new Unit instance in the given system.
   */
  toSystem(system) {
    // Go through the currently defined units, identify their dimension and look
    // for the corresponding dimension in the given unit system. If one is
    // present, convert to it. Otherwise convert to base dimensions.
    const UNITS = UnitMathJS.UNITS
    const PREFIXES = UnitMathJS.PREFIXES
    const BASE_DIMENSIONS = UnitMathJS.BASE_DIMENSIONS
    const BASE_UNITS = UnitMathJS.BASE_UNITS
    const proposedUnitList = []
    for (const unit of this.mathjsUnit.units) {
      const dimension = unit.unit.base.key
      const newUnitName = system?.[dimension]?.name
      const newUnit = UNITS[newUnitName]
      // If the unit for this dimension is defined, use it
      if (!isNil(newUnitName)) {
        proposedUnitList.push({
          unit: newUnit,
          prefix: PREFIXES.NONE[''],
          power: unit.power || 0
        })
      // Otherwise convert to base units
      } else {
        let missingBaseDim = false
        const baseUnit = BASE_UNITS[dimension]
        const newDimensions = baseUnit.dimensions
        for (let i = 0; i < BASE_DIMENSIONS.length; i++) {
          const baseDim = BASE_DIMENSIONS[i]
          if (Math.abs(newDimensions[i] || 0) > 1e-12) {
            if (has(system, baseDim)) {
              proposedUnitList.push({
                unit: UNITS[system[baseDim].name],
                prefix: PREFIXES.NONE[''],
                power: unit.power ? newDimensions[i] * unit.power : 0
              })
            } else {
              missingBaseDim = true
            }
          }
        }
        if (missingBaseDim) {
          throw Error(`The given unit system does not contain the required unit definitions for converting ${unit.name} with dimension ${dimension}.`)
        }
      }
    }

    const ret = this.mathjsUnit.clone()
    ret.units = proposedUnitList
    return new Unit(ret)
  }
}

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
    return deepMap(value, (x) => this.unit.mathjsUnit._normalize(x))
  }

  /**
   * Convert value from base units to currently set units.
   * @param {n-dimensional array} value Value in base units.
   * @returns Value in currently set units.
   */
  denormalize(value) {
    return deepMap(value, (x) => this.unit.mathjsUnit._denormalize(x))
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
}

/**
 * Convenience function for getting compatible units for a given dimension.
 * Returns all compatible units that have been registered.
 *
 * @param {string} dimension The dimension.
 * @returns Array of compatible units.
 */
export function getUnits(dimension) {
  return dimensionMap?.[dimension]?.units || []
}

/**
 * Convenience function for parsing unit information from a string.
 *
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
    return {value, valueString, unitString, error: 'Enter a valid unit'}
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
    return {valueString, value, unitString, error: `Unknown unit '${unitString}'`}
  }
  if (unit.dimension(false) !== dimension) {
    return {valueString, value, unitString, unit, error: `Incompatible unit`}
  }
  return {value, valueString, unit, unitString}
}
