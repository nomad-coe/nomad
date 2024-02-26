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

import React, {useState, useMemo, useContext, useCallback, useEffect} from 'react'
import PropTypes from 'prop-types'
import { isNil, isFunction, startCase, toLower, cloneDeep } from 'lodash'
import { Unit as UnitMathJS, createUnit } from 'mathjs'
import { unitList, unitPrefixes as prefixes } from '../../config'

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
export const unitToAbbreviationMap = {}
const unitDefinitions = {}
for (let def of unitList) {
  const name = def.name
  def = {
    ...def,
    baseName: def.dimension,
    prefixes: 'pint'
  }
  unitDefinitions[name] = def
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

// Export unit options for each unit and dimension
export const unitMap = Object.fromEntries(unitList.map(x => [x.name, x]))
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
 * React context for interacting with unit configurations.
 */
export const unitContext = React.createContext()
export const UnitProvider = React.memo(({initialUnitSystems, initialSelected, children}) => {
  const resetUnitSystems = useState(cloneDeep(initialUnitSystems))[0]
  const [isReset, setIsReset] = useState(true)
  const [unitSystems, setUnitSystems] = useState(cloneDeep(initialUnitSystems))
  const [selectedSystem, setSelectedSystem] = useState(initialSelected)
  const [currentScope, setCurrentScope] = useState()
  const [scopes, setScopes] = useState({})
  const [scope, setScope] = useState()

  const reset = useCallback(() => {
    setUnitSystems(cloneDeep(resetUnitSystems))
    setIsReset(true)
    if (currentScope && scopes?.[currentScope]?.system) {
      setSelectedSystem(scopes[currentScope].system)
    }
  }, [currentScope, scopes, resetUnitSystems])

  const setScopeInfo = useCallback((scopeKey, unit) => {
    if (!(scopeKey in scopes)) {
      scopes[scopeKey] = {system: unit || initialSelected, unitSystems: unitSystems}
      setScopes({...scopes})
    }
    setCurrentScope(scopeKey)
  }, [scopes, initialSelected, unitSystems])

  useEffect(() => {
    if (currentScope in scopes) {
      setSelectedSystem(scopes[currentScope].system)
      setUnitSystems(scopes[currentScope].unitSystems)
    }
  }, [currentScope, scopes])

  const setSelected = useCallback((unitSystem) => {
    setIsReset(false)
    setSelectedSystem(unitSystem)
    setScopes(scopes => {
      scopes[currentScope].system = unitSystem
      return {...scopes}
    })
  }, [currentScope])

  const setUnits = useCallback((value) => {
    setIsReset(false)
    setUnitSystems(old => {
      const newSystems = {...old}
      newSystems[selectedSystem].units = isFunction(value)
        ? value(newSystems[selectedSystem].units)
        : value
      setScopes(scopes => {
        scopes[currentScope].unitSystems = newSystems
        return {...scopes}
      })
      return newSystems
    })
  }, [currentScope, selectedSystem])

  const values = useMemo(() => {
    return {
      units: unitSystems[selectedSystem].units,
      setUnits: setUnits,
      unitSystems,
      unitMap,
      dimensionMap,
      selected: selectedSystem,
      setSelected,
      reset,
      isReset,
      setScopeInfo,
      scope,
      setScope
    }
  }, [unitSystems, selectedSystem, setUnits, setSelected, reset, isReset, setScopeInfo, scope, setScope])

  return <unitContext.Provider value={values}>
    {children}
  </unitContext.Provider>
})

UnitProvider.propTypes = {
  initialUnitSystems: PropTypes.object,
  initialSelected: PropTypes.string,
  children: PropTypes.node
}

/**
 * Convenience hook for using the current unit context.
 *
 * @returns Object containing the currently set units for each dimension (e.g.
 * {energy: 'joule'})
 */
export const useUnitContext = (scopeKey, defaultUnitSystem) => {
  const context = useContext(unitContext)
  const {setScopeInfo, scope, setScope, ...others} = context

  useEffect(() => {
    setScope(!scopeKey ? scopeKey : 'Schema')
  }, [scopeKey, setScope])

  useEffect(() => {
    if (scopeKey && scope === 'Schema') {
      setScopeInfo(scopeKey, defaultUnitSystem)
    } else {
      setScopeInfo('Global')
    }
  }, [scopeKey, scope, setScopeInfo, defaultUnitSystem])

  return {...others, scope, setScope}
}
