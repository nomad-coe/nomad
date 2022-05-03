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
import React, {useCallback, useEffect, useMemo, useRef, useState} from 'react'
import {TextField, makeStyles, Box, MenuItem, Checkbox, Tooltip} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import {debounce} from 'lodash'
import {TextFieldWithHelp, getFieldProps} from './StringEditQuantity'
import {useErrors} from '../errors'

export const NumberField = React.memo((props) => {
  const {onChange, value, dataType, minValue, maxValue, displayUnit, unit, convertInPlace, ...otherProps} = props
  const inputRef = useRef()
  const dimension = unit && unitMap[unit].dimension
  const units = dimension && conversionMap[dimension]?.units
  const requestedUnit = useRef(displayUnit) // Changing variables from inside useCallback will be lost after each render.

  const getValue = useCallback((value) => {
    if (value === undefined) return ''
    if (units === undefined) {
      requestedUnit.current = displayUnit
      return value
    }
    const numberPart = value.match(/^[+-]?((\d+\.\d+|\d+\.|\.\d?|\d+)(e|e\+|e-)\d+|(\d+\.\d+|\d+\.|\.\d?|\d+))?/)
    if (numberPart !== undefined) {
      const unitPart = value.substring(numberPart[0].length).trim().toLowerCase()
      let unitObject
      try {
        unitObject = new Unit(unitPart)
      } catch {}
      if (unitObject !== undefined && units.find(value => unitObject.unitDef.name === value)) {
        requestedUnit.current = unitObject.unitDef.name
        return numberPart[0]
      }
    }
    requestedUnit.current = displayUnit
    return value
  }, [displayUnit, units])

  const createInputValue = useCallback(value => {
    if (value === undefined) return ''
    if (isNaN(value)) return ''

    if (unit && requestedUnit.current) {
      value = convertUnit(value, unit, requestedUnit.current)
    }
    // Make sure that the formatting is not overwriting the format used to enter the value
    const inputValue = inputRef.current?.value || ''
    if (value === Number(getValue(inputValue))) {
      return value
    }
    return String(value) // TODO when to use toLocaleString, when to use exp. format?
  }, [unit, requestedUnit, getValue])

  const [inputValue, setInputValue] = useState(createInputValue(value))
  const [error, setError] = useState('')

  useEffect(() => {
    if (convertInPlace) {
      setInputValue(createInputValue(value))
    }
  }, [createInputValue, convertInPlace, value, getValue])

  const checkAndGetNumberValue = useCallback((value) => {
    value = getValue(value)
    if (['int64', 'int32', 'int'].includes(dataType)) {
      const num = Number(value)
      return [Number.isInteger(num), num]
    } else if (['uint64', 'uint32', 'uint'].includes(dataType)) {
      const num = Number(value)
      return [Number.isInteger(num) && num >= 0, num]
    } else if (['float64', 'float32', 'float'].includes(dataType)) {
      const num = Number(value)
      return [!isNaN(num), num]
    }
  }, [dataType, getValue])

  const handleChange = useCallback(() => {
    if (!inputRef.current) return
    const value = inputRef.current.value.trim().replace(/,/g, '.')
    if (value === '') {
      if (onChange) {
        onChange(undefined, requestedUnit.current)
      }
      setError('')
      return
    }

    let [isNumber, number] = checkAndGetNumberValue(value)
    if (!isNumber) {
      setError('Enter a valid value.')
      return
    }

    if (!isNaN(maxValue) && maxValue < number) {
      setError(`Enter a value smaller than ${maxValue}.`)
      return
    }

    if (!isNaN(minValue) && minValue > number) {
      setError(`Enter a value larger than ${minValue}.`)
      return
    }

    if (requestedUnit.current && unit) {
      number = convertUnit(number, requestedUnit.current, unit)
    }

    if (!convertInPlace) {
      setInputValue(getValue(value))
    }

    if (onChange) {
      onChange(number, requestedUnit.current)
    }
    setError('')
  }, [checkAndGetNumberValue, maxValue, minValue, unit, convertInPlace, onChange, getValue])

  const debouncedHandleChange = useMemo(() => {
    return debounce(handleChange, 500)
  }, [handleChange])

  const handleBlur = useCallback(() => {
    debouncedHandleChange.flush()
  }, [debouncedHandleChange])

  const handleInputChange = useCallback(event => {
    setInputValue(event.target.value)
    debouncedHandleChange()
  }, [debouncedHandleChange])

  return <TextFieldWithHelp
    ref={inputRef}
    fullWidth variant='filled' size='small'
    value={String(inputValue)}
    onBlur={handleBlur} error={!!error} helperText={error}
    onChange={handleInputChange}
    data-testid='number-edit-quantity-value'
    {...otherProps}
  />
})
NumberField.propTypes = {
  onChange: PropTypes.func.isRequired,
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  dataType: PropTypes.string,
  value: PropTypes.number,
  displayUnit: PropTypes.string,
  unit: PropTypes.string,
  convertInPlace: PropTypes.bool
}

export const NumberEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, defaultDisplayUnit, ...otherProps} = props
  const systemUnits = useUnits()
  const {raiseError} = useErrors()
  const hasUnit = quantityDef.unit
  const dimension = hasUnit && unitMap[quantityDef.unit].dimension
  const [checked, setChecked] = useState(true)
  let defaultUnit
  if (defaultDisplayUnit) {
    try {
      defaultUnit = new Unit(defaultDisplayUnit)
      if (!conversionMap[dimension].units.find(value => value === defaultUnit.unitDef.name)) {
        raiseError(`The provided defaultDisplayUnit for ${quantityDef.name} is not related to this field.`)
      }
    } catch (e) {
      raiseError(`The provided defaultDisplayUnit for ${quantityDef.name} field is not valid.`)
    }
  }
  const [unit, setUnit] = useState((defaultUnit && defaultUnit.unitDef.name) || systemUnits[dimension] || quantityDef.unit)

  const handleChange = useCallback((value, unit) => {
    onChange(value)
    setUnit(unit)
  }, [onChange])

  const handleUnitChange = useCallback((newUnit) => {
    if (!checked && quantityDef.unit && newUnit) {
      let displayedValue = convertUnit(value, quantityDef.unit, unit)
      let storedValue = convertUnit(displayedValue, newUnit, quantityDef.unit)
      onChange(storedValue)
    }
    setUnit(newUnit)
  }, [checked, onChange, quantityDef.unit, unit, value])

  return <Box display='flex'>
    <NumberField
      value={value} onChange={handleChange}
      unit={quantityDef.unit}
      convertInPlace={checked}
      displayUnit={unit}
      dataType={quantityDef.type?.type_data}
      {...getFieldProps(quantityDef)}
      {...otherProps}
    />
    {hasUnit && (
      <Box display='flex'>
        <Tooltip title={'If checked, numeric value is converted when the unit is changed.'}>
          <Checkbox
            onChange={event => setChecked(event.target.checked)}
            color="primary"
            checked={checked}
          />
        </Tooltip>
        <UnitSelect defaultUnit={quantityDef.unit} unit={unit} onChange={handleUnitChange}/>
      </Box>
    )}
  </Box>
})
NumberEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.number,
  onChange: PropTypes.func,
  defaultDisplayUnit: PropTypes.string
}

export const useUnitSelectStyles = makeStyles(theme => ({
  root: {
    width: '150px'
  }
}))

export const UnitSelect = React.memo(({defaultUnit, unit, onChange}) => {
  const classes = useUnitSelectStyles()
  const dimension = unitMap[defaultUnit].dimension
  const units = conversionMap[dimension].units

  const handleUnitChange = useCallback(event => {
    onChange(event.target.value)
  }, [onChange])

  return (
    <TextField
      className={classes.root} variant='filled' size='small' select
      label="unit" value={unit}
      onChange={handleUnitChange}
      data-testid = 'number-edit-quantity-unit'
    >
      {units.map(unit => (
        <MenuItem key={unit} value={unit}>
          {(new Unit(unit)).label}
        </MenuItem>
      ))}
    </TextField>
  )
})
UnitSelect.propTypes = {
  defaultUnit: PropTypes.string.isRequired,
  unit: PropTypes.string.isRequired,
  onChange: PropTypes.func
}
