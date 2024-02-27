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
import React, { useState, useEffect, useCallback, useRef } from 'react'
import PropTypes from 'prop-types'
import {useUnitContext} from './UnitContext'
import UnitInput from './UnitInput'

/**
 * Controls the unit for the specified dimension in the current unit context.
 */
const UnitDimensionSelect = React.memo(({label, dimension, onChange, disabled}) => {
  const {units, setUnits} = useUnitContext()
  const [error, setError] = useState()
  const unit = units?.[dimension]
  const disabledFinal = disabled || unit?.locked
  const labelFinal = label || dimension
  const oldValue = useRef(unit?.definition)
  const [inputValue, setInputValue] = useState(unit?.definition)

  // React to changes in units
  useEffect(() => {
    setError(undefined)
    setInputValue(unit?.definition)
  }, [unit, dimension])

  const handleAccept = useCallback((unitString, unit) => {
    setUnits(old => {
      const newUnits = {
        ...old,
        [dimension]: {...old[dimension], definition: unitString}
      }
      return newUnits
    })
    onChange?.(unit)
    oldValue.current = unitString
  }, [dimension, onChange, setUnits])

  const handleSelect = useCallback((unitString, unit) => {
    handleAccept(unit.label(), unit)
  }, [handleAccept])

  const handleChange = useCallback((value) => {
    oldValue.current = value
    setInputValue(value)
  }, [])

  return (unit && dimension !== 'dimensionless')
    ? <UnitInput
      value={inputValue}
      label={labelFinal}
      onChange={handleChange}
      onAccept={handleAccept}
      onSelect={handleSelect}
      onError={setError}
      dimension={dimension}
      error={error}
      disabled={disabledFinal}
      disableGroup
    />
    : null
})
UnitDimensionSelect.propTypes = {
  value: PropTypes.string,
  label: PropTypes.string,
  dimension: PropTypes.string,
  onChange: PropTypes.func,
  disabled: PropTypes.bool
}

export default UnitDimensionSelect
