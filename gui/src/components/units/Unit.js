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
import {isNil, has, isString} from 'lodash'
import {Unit as UnitMathJS} from 'mathjs'
import {unitToAbbreviationMap} from './UnitContext'

export const DIMENSIONLESS = 'dimensionless'

/**
 * Helper class for persisting unit information.
 *
 * Builds upon the math.js Unit class, but adds additional functionality,
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
      unit = normalizeExpression(unit)
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
   * Checks if the given unit has the same base dimensions as this one.
   * @param {str | Unit} unit Unit to compare to
   * @returns boolean Whether the units have the same base dimensions.
   */
  equalBase(unit) {
    if (isString(unit)) {
      unit = normalizeExpression(unit)
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
      if (unit.base.key === DIMENSIONLESS) return ''
      return abbreviate
        ? unitToAbbreviationMap?.[unit.name] || unit.name
        : unit.name
    }

    function getPrefix(unit, original) {
      if (!abbreviate) return original
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
   * @param {boolean} base Whether to return dimension in base units. Otherwise
   * the original unit dimensions are used.
   * @returns The dimensionality as a string, e.g. 'time^2 energy mass^-2'
   */
  dimension(base = false) {
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
    const dims = Object.entries(dimensionMap).filter(d => d[1] !== 0)
    return dims.length > 0
      ? dims.map(d => `${d[0]}${((d[1] < 0 || d[1] > 1) && `^${d[1]}`) || ''}`).join(' ')
      : DIMENSIONLESS
  }

  /**
   * Function for converting to another unit.
   *
   * @param {str | Unit} unit The target unit
   * @returns A new Unit expressed in the given units.
   */
  to(unit) {
    if (isString(unit)) {
      unit = normalizeExpression(unit)
    } else if (unit instanceof Unit) {
      unit = unit.label()
    } else {
      throw Error('Unknown unit type. Please provide the unit as as string or as instance of Unit.')
    }

    // We cannot directly feed the unit string into Math.js, because it will try
    // to parse units like 1/<unit> as Math.js units which have values, and then
    // will raise an exception when converting between valueless and valued
    // unit. The workaround is to explicitly define a valueless unit.
    unit = new UnitMathJS(undefined, unit === '' ? DIMENSIONLESS : unit)
    return new Unit(this.mathjsUnit.to(unit))
  }

  /**
   * Function for converting the value of this Unit to the SI unit system.
   *
   * @returns A new Unit instance in the SI unit system.
   */
  toSI() {
    return this.toSystem({
      "dimensionless": { "definition": "dimensionless" },
      "length": { "definition": "m" },
      "mass": { "definition": "kg" },
      "time": { "definition": "s" },
      "current": { "definition": "A" },
      "temperature": { "definition": "K" },
      "luminosity": { "definition": "cd" },
      "luminous_flux": { "definition": "lm" },
      "substance": { "definition": "mol" },
      "angle": { "definition": "rad" },
      "information": { "definition": "bit" },
      "force": { "definition": "N" },
      "energy": { "definition": "J" },
      "power": { "definition": "W" },
      "pressure": { "definition": "Pa" },
      "charge": { "definition": "C" },
      "resistance": { "definition": "Ω" },
      "conductance": { "definition": "S" },
      "inductance": { "definition": "H" },
      "magnetic_flux": { "definition": "Wb" },
      "magnetic_field": { "definition": "T" },
      "frequency": { "definition": "Hz" },
      "luminance": { "definition": "nit" },
      "illuminance": { "definition": "lx" },
      "electric_potential": { "definition": "V" },
      "capacitance": { "definition": "F" },
      "activity": { "definition": "kat" }
    })
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
      const newUnitDefinition = system?.[dimension]?.definition
      // If the unit for this dimension is defined, use it
      if (!isNil(newUnitDefinition)) {
        const newUnit = new Unit(newUnitDefinition)
        for (const unitDef of newUnit.mathjsUnit.units) {
          proposedUnitList.push({...unitDef, power: unitDef.power * unit.power})
        }
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
                unit: UNITS[system[baseDim].definition],
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
 * Normalizes the given expression into a format that can be parsed by MathJS.
 *
 * This function will replace the Pint power symbol of '**' with the symbol
 * '^' used by MathJS. In addition, we convert any 'delta'-units (see:
 * https://pint.readthedocs.io/en/stable/nonmult.html) into their regular
 * counterparts: MathJS will automatically ignore the offset when using
 * non-multiplicative units in expressions.
 *
 * @param {str} expression Expression
 * @returns string Expression in normalized form
 */
export function normalizeExpression(expression) {
  let normalized = expression.replace(/\*\*/g, '^')
  normalized = normalized.replace(/delta_/g, '')
  normalized = normalized.replace(/Δ/g, '')
  return normalized
}
