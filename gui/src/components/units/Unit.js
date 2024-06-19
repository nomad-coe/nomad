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
import {parseInternal} from './common'

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
      unit = parseInternal(unit, {requireUnit: true}).unit
    } else if (unit instanceof Unit) {
      unit = unit.mathjsUnit.clone()
    } else if (unit instanceof UnitMathJS) {
      unit = unit.clone()
    } else {
      throw Error('Please provide the unit as a string or as an instance of Unit.')
    }
    this.mathjsUnit = unit
  }

  /**
   * Checks if the given unit has the same base dimensions as this one.
   * @param {str | Unit} unit Unit to compare to
   * @returns boolean Whether the units have the same base dimensions.
   */
  equalBase(unit) {
    if (isString(unit)) {
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
  label(abbreviate = true, showDelta = false) {
    // TODO: The label caching is disabled for now. Because Quantities are
    // stored as recoil.js atoms, they become immutable which causes problems
    // with internal state mutation.
    const units = this.mathjsUnit.units
    let strNum = ''
    let strDen = ''
    let nNum = 0
    let nDen = 0

    function getDelta(unit) {
      return (showDelta && unit.delta)
        ? abbreviate ? 'Δ' : 'delta_'
        : ''
    }

    function getName(unit) {
      if (unit.base.key === DIMENSIONLESS) return ''
      return abbreviate
        ? unitToAbbreviationMap?.[unit.name] || unit.name
        : unit.name
    }

    function getPrefix(original) {
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
      const unitDef = units[i]
      if (unitDef.power > 0) {
        nNum++
        const prefix = getPrefix(unitDef.prefix.name)
        const name = getName(unitDef.unit)
        const delta = getDelta(unitDef)
        strNum += ` ${delta}${prefix}${name}`
        if (Math.abs(unitDef.power - 1.0) > 1e-15) {
          strNum += '^' + unitDef.power
        }
      } else if (unitDef.power < 0) {
        nDen++
      }
    }

    if (nDen > 0) {
      for (let i = 0; i < units.length; i++) {
        const unitDef = units[i]
        if (unitDef.power < 0) {
          const prefix = getPrefix(unitDef.prefix.name)
          const name = getName(unitDef.unit)
          const delta = getDelta(unitDef)
          if (nNum > 0) {
            strDen += ` ${delta}${prefix}${name}`
            if (Math.abs(unitDef.power + 1.0) > 1e-15) {
              strDen += '^' + (-unitDef.power)
            }
          } else {
            strDen += ` ${delta}${prefix}${name}`
            strDen += '^' + (unitDef.power)
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
    let mathjsUnit
    if (unit instanceof Unit) {
      mathjsUnit = unit.mathjsUnit
    } else if (isString(unit)) {
      mathjsUnit = parseInternal(unit, {requireUnit: true}).unit
    } else {
      throw Error('Unknown unit type. Please provide the unit as as string or as instance of Unit.')
    }
    return new Unit(this.mathjsUnit.to(mathjsUnit))
  }

  /**
   * Converts all units to their delta versions.
   *
   * @returns A new Unit with delta units.
   */
  toDelta() {
    const unitMathJS = this.mathjsUnit.clone()
    for (const unit of unitMathJS.units) {
      unit.delta = true
    }
    return new Unit(unitMathJS)
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
   * If the converted unit has any delta-units, the converted units will also
   * become delta-units.
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
          proposedUnitList.push({...unitDef, power: unitDef.power * unit.power, delta: unit.delta})
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
              let unitDef
                const baseUnitDefinition = system[baseDim].definition
                // First try to look up the unit in the list of known units:
                // most likely it is already there. If not found, then we need
                // to try and parse the unit definition: this happens e.g. when
                // units with prefixes are used as base units in the unit
                // system.
                unitDef = UNITS[system[baseDim]?.definition]
                if (!unitDef) {
                  try {
                    const baseUnit = new Unit(baseUnitDefinition)
                    const baseUnitList = baseUnit.mathjsUnit.units
                    if (baseUnitList.length === 1) unitDef = baseUnitList[0].unit
                  } catch (e) {}
                }
              if (!unitDef) {
                throw Error(`Could not find base unit for "${system[baseDim].definition}"`)
              }
              proposedUnitList.push({
                unit: unitDef,
                prefix: PREFIXES.NONE[''],
                power: unit.power ? newDimensions[i] * unit.power : 0,
                delta: unit.delta
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
