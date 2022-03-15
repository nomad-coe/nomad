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
import {TextField, makeStyles, Box, MenuItem} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import {debounce} from 'lodash'
import {TextFieldWithHelp, getFieldProps} from './StringEditQuantity'

export const NumberField = React.memo((props) => {
  const {onChange, value, dataType, minValue, maxValue, displayUnit, unit, ...otherProps} = props
  const inputRef = useRef()

  const createInputValue = useCallback(value => {
    if (isNaN(value)) {
      return ''
    }
    if (unit && displayUnit) {
      value = convertUnit(value, unit, displayUnit)
    }
    // Make sure that the formatting is not overwriting the format used to enter the value
    const inputValue = inputRef.current?.value || ''
    if (value === Number(inputValue)) {
      return inputValue
    }
    return String(value) // TODO when to use toLocaleString, when to use exp. format?
  }, [unit, displayUnit])

  const [inputValue, setInputValue] = useState(createInputValue(value))
  const [error, setError] = useState('')

  useEffect(() => {
    setInputValue(createInputValue(value))
  }, [createInputValue, value])

  const checkAndGetNumberValue = useCallback((value) => {
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
  }, [dataType])

  const handleChange = useCallback(() => {
    const value = inputRef.current.value.trim().replace(/,/g, '.')
    if (value === '') {
      if (onChange) {
        onChange(undefined)
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

    if (displayUnit && unit) {
      number = convertUnit(number, displayUnit, unit)
    }
    if (onChange) {
      onChange(number)
    }
    setError('')
  }, [maxValue, minValue, onChange, checkAndGetNumberValue, displayUnit, unit, inputRef])

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
    value={inputValue}
    onBlur={handleBlur} error={!!error} helperText={error}
    onChange={handleInputChange}
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
  unit: PropTypes.string
}

export const NumberEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, ...otherProps} = props
  const systemUnits = useUnits()
  const hasUnit = quantityDef.unit
  const dimension = hasUnit && unitMap[quantityDef.unit].dimension
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)

  return <Box display='flex'>
    <NumberField
      value={value} onChange={onChange}
      unit={quantityDef.unit}
      displayUnit={unit}
      dataType={quantityDef.type?.type_data}
      {...getFieldProps(quantityDef)}
      {...otherProps}
    />
    {hasUnit && (
      <UnitSelect defaultUnit={quantityDef.unit} unit={unit} onChange={setUnit}/>
    )}
  </Box>
})
NumberEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.number,
  onChange: PropTypes.func
}

export const useUnitSelectStyles = makeStyles(theme => ({
  root: {
    marginLeft: theme.spacing(1),
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
