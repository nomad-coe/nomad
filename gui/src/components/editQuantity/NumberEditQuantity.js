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
import {TextField, makeStyles, Box, Checkbox, Tooltip} from '@material-ui/core'
import Autocomplete from '@material-ui/lab/Autocomplete'
import PropTypes from 'prop-types'
import {Quantity, parseQuantity} from '../units/Quantity'
import {Unit} from '../units/Unit'
import {getUnits} from '../units/UnitContext'
import {debounce, isNil} from 'lodash'
import {TextFieldWithHelp, getFieldProps} from './StringEditQuantity'
import {useDisplayUnit} from "../units/useDisplayUnit"
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const NumberField = React.memo((props) => {
  const {onChange, onInputChange, dimension, value, dataType, minValue, unit, maxValue, displayUnit, convertInPlace, debounceTime, ...otherProps} = props
  const [error, setError] = useState('')
  const previousValue = useRef('')
  const previousNumber = useRef()
  const previousUnitLabel = useRef(displayUnit?.label())
  const previousNumberPart = useRef('')
  const [inputValue, setInputValue] = useState('')

  const fixDigits = useCallback((value) => {
    return value.toPrecision(15)
      .replace(/0+(e|$)/, '$1')
      .replace(/\.$/, '')
      .replace(/\.e/, 'e')
  }, [])

  const getStoredValue = useCallback((value, displayUnit) => {
    const valueInBaseUnit = (isNil(value) || isNil(unit))
      ? value
      : new Quantity(value, displayUnit).to(unit).value()
    return valueInBaseUnit && fixDigits(Number(valueInBaseUnit))
  }, [fixDigits, unit])

  // Whenever a new value arrives or the units change, change the text
  // accordingly. If the field is associated with a unit and the same numerical
  // value as set previously comes back with the same unit, preserve the
  // formatting by reusing the saved input text. Fields that do not have a unit
  // simply take in the new value.
  useEffect(() => {
    if (!isNil(unit)) {
      if (convertInPlace) {
        let newVal
        if (!isNil(value) && value === previousNumber.current && displayUnit?.label() === previousUnitLabel.current) {
          newVal = previousNumberPart.current
        } else if (isNil(value) || isNaN(value)) {
          newVal = ''
        } else if (unit && displayUnit) {
          newVal = new Quantity(value, unit).to(displayUnit).value()
        }
        previousUnitLabel.current = displayUnit?.label()
        previousNumber.current = value
        previousNumberPart.current = newVal
        const inputValue = newVal && fixDigits(Number(newVal))
        previousValue.current = inputValue
        setInputValue(inputValue)
        onInputChange(inputValue)
      }
    } else {
      const inputValue = (isNil(value) || isNaN(value)) ? '' : String(value)
      previousValue.current = inputValue
      setInputValue(inputValue)
    }
  }, [value, displayUnit, unit, convertInPlace, onInputChange, fixDigits])

  // Given the text input, returns a the extracted number, unit and possible
  // errors.
  const parseInput = useCallback((input) => {
    if (isNil(input) || input === '') {
      previousNumberPart.current = ''
      return {}
    }

    // Try to parse the quantity. Value is required, unit is optional.
    const {unit: parsedUnit, value, valueString, error} = parseQuantity(input, dimension, true, false)
    previousNumberPart.current = valueString
    if (parsedUnit) {
      previousUnitLabel.current = parsedUnit.label()
    }
    if (error) {
      return {error}
    }

    // Validate number datatype
    let isValid = true
    if (['int64', 'int32', 'int'].includes(dataType)) {
      isValid = Number.isInteger(value)
    } else if (['uint64', 'uint32', 'uint'].includes(dataType)) {
      isValid = Number.isInteger(value) && value >= 0
    } else if (['float64', 'float32', 'float'].includes(dataType)) {
      isValid = !isNaN(value)
    }
    if (!isValid) {
      return {error: `Enter a valid numerical value with datatype '${dataType}'`}
    }

    // Validate number range
    const storingValue = Number(getStoredValue(value, parsedUnit || displayUnit))
    if (!isNaN(storingValue)) {
      if (!isNaN(maxValue) && maxValue < storingValue) {
        return {error: `Enter a value that is equal or smaller than ${maxValue}` + (unit ? ` (${unit})` : '')}
      }
      if (!isNaN(minValue) && minValue > storingValue) {
        return {error: `Enter a value that is equal or larger than ${minValue}` + (unit ? ` (${unit})` : '')}
      }
    }

    return {value: value, unit: parsedUnit}
  }, [dataType, dimension, displayUnit, getStoredValue, maxValue, minValue, unit])

  // The final handler for text input change
  const handleChange = useCallback((inputValue) => {
    // If the input has not been changed meaningfully, do nothing.
    const trimmed = inputValue?.trim()
    if (trimmed === previousValue.current) {
      return
    }

    // Store the new value in a ref, reset errors
    previousValue.current = trimmed
    setError('')

    // Try to interpret the value, if a meaningful value was returned, send the
    // new value in a normalized form to the handler.
    const parsedInput = parseInput(trimmed)
    const {value: number, error} = parsedInput
    const newUnit = parsedInput.unit || displayUnit

    if (error) {
      setError(error)
    } else {
      if (!convertInPlace) {
        setInputValue(previousNumberPart.current)
        onInputChange(previousNumberPart.current)
      }
      if (onChange) {
        const storedValue = getStoredValue(number, newUnit)
        previousNumber.current = storedValue
        const numberedStoredValue = Number(storedValue)
        onChange(!isNaN(numberedStoredValue) ? numberedStoredValue : undefined, newUnit)
        onInputChange(number)
      }
    }
  }, [parseInput, displayUnit, convertInPlace, onChange, onInputChange, getStoredValue])

  // Routes text field changes to a handler after a debounce time
  const debouncedHandleChange = useMemo(() => debounce(handleChange, 500), [handleChange])

  // When input changes, saves the entered text into a state and queues a
  // processing event with debounce time.
  const handleInputChange = useCallback(event => {
    setInputValue(event.target.value)
    if (debounceTime) {
      debouncedHandleChange(event.target.value)
    }
  }, [debounceTime, debouncedHandleChange])

  // Request immediate validation on blur
  const handleBlur = useCallback(() => {
    handleChange(inputValue)
  }, [handleChange, inputValue])

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    value={String(inputValue)}
    onBlur={handleBlur}
    error={!!error}
    helperText={error}
    onChange={handleInputChange}
    onKeyDown={(event) => { if (event.key === 'Enter') handleChange(inputValue) }}
    data-testid='number-edit-quantity-value'
    {...otherProps}
  />
})
NumberField.propTypes = {
  onChange: PropTypes.func.isRequired,
  onInputChange: PropTypes.func.isRequired,
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  dataType: PropTypes.string,
  value: PropTypes.number,
  unit: PropTypes.string,
  displayUnit: PropTypes.object,
  dimension: PropTypes.string,
  convertInPlace: PropTypes.bool,
  debounceTime: PropTypes.number
}

export const NumberEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, ...otherProps} = props
  const defaultUnit = useMemo(() => quantityDef.unit && new Unit(quantityDef.unit), [quantityDef])
  const dimension = defaultUnit && defaultUnit.dimension(false)
  const dimensionBase = defaultUnit && defaultUnit.dimension(true)
  const [checked, setChecked] = useState(true)
  const [displayedValue, setDisplayedValue] = useState(true)
  const {defaultDisplayUnit: deprecatedDefaultDisplayUnit, ...fieldProps} = getFieldProps(quantityDef)
  const displayUnit = useDisplayUnit(quantityDef)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)
  const [unit, setUnit] = useState(displayUnit)

  // Set the unit if display unit changes
  useEffect(() => {
    setUnit(displayUnit)
  }, [displayUnit])

  // Get a list of unit options for this field
  const options = useMemo(() => {
    const units = [...getUnits(dimension)].map(x => new Unit(x).label())
    unit && units.push(unit.label())
    defaultUnit && units.push(defaultUnit.label())
    displayUnit && units.push(displayUnit.label())
    return [...new Set(units)]
  }, [displayUnit, defaultUnit, dimension, unit])

  // Handle a change in NumberField input
  const handleChange = useCallback((value, unit) => {
    onChange(value)
    setUnit(unit)
  }, [onChange])

  // Handle a change in NumberField input
  const handleInputChange = useCallback((displayedValue) => {
    setDisplayedValue(displayedValue)
  }, [])

  // Handle a change in the unit dialog
  const handleUnitChange = useCallback((newUnit) => {
    if (!checked && quantityDef.unit && newUnit && !isNil(value)) {
      const storedValue = new Quantity(Number(displayedValue), newUnit).to(quantityDef.unit).value()
      onChange(storedValue)
    }
    setUnit(new Unit(newUnit))
  }, [checked, displayedValue, onChange, quantityDef.unit, value])

  return <Box display='flex'>
    <NumberField
      value={value}
      onChange={handleChange}
      onInputChange={handleInputChange}
      dimension={dimension}
      convertInPlace={checked}
      unit={quantityDef?.unit}
      displayUnit={unit}
      dataType={quantityDef.type?.type_data}
      {...fieldProps}
      {...otherProps}
      label={label}
    />
    {unit && (
      <Box display='flex'>
        <Tooltip title={'If checked, numeric value is converted when the unit is changed.'}>
          <Checkbox
            onChange={event => setChecked(event.target.checked)}
            color="primary"
            checked={checked}
          />
        </Tooltip>
        <UnitSelect options={options} unit={unit} dimension={dimensionBase} onChange={handleUnitChange}/>
      </Box>
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
    width: '150px'
  }
}))

export const UnitSelect = React.memo(({options, unit, onChange, dimension, disabled}) => {
  const classes = useUnitSelectStyles()
  const [error, setError] = useState()
  const [value, setValue] = useState(unit.label())
  const [open, setOpen] = useState(false)

  // Set the input value when unit changes
  useEffect(() => {
    unit && setValue(unit.label())
  }, [unit])

  // Validate input and submit unit if valid
  const submit = useCallback((val) => {
    const {unit, error} = parseQuantity(val, dimension, false, true)
    if (error) {
      setError(error)
    } else {
      onChange(unit)
    }
  }, [dimension, onChange])

  const handleInputChange = useCallback((event, value) => {
    setError()
    setValue(value)
  }, [])

  const handleOptionChange = useCallback((event, value) => {
    setError()
    submit(value)
  }, [submit])

  const finalOptions = useMemo(() => {
    return options
  }, [options])

  // Submit on enter
  const handleEnter = useCallback((event) => {
    if (event.key === 'Enter') {
      setOpen(false)
      submit(value)
      event.stopPropagation()
      event.preventDefault()
    }
  }, [submit, value])

  return <Autocomplete
    className={classes.root}
    size='small'
    freeSolo
    disabled={disabled}
    forcePopupIcon
    disableClearable
    value={value}
    options={finalOptions}
    getOptionLabel={(option) => option}
    filterOptions={(options) => options}
    onInputChange={handleInputChange}
    onChange={handleOptionChange}
    open={open}
    onOpen={() => { if (value.trim() !== '') { setOpen(true) } }}
    onClose={() => setOpen(false)}
    onBlur={() => submit(value)}
    renderInput={(params) => <TextField
      {...params}
      label="Unit"
      data-testid = 'number-edit-quantity-unit'
      variant="filled"
      error={!!error}
      helperText={error}
      onKeyDown={handleEnter}
      InputProps={{
        ...params.InputProps
        // endAdornment: <HelpAdornment title={'Moi'} description={'Terve'}/>
      }}
    />}
  />
})
UnitSelect.propTypes = {
  options: PropTypes.arrayOf(PropTypes.string),
  unit: PropTypes.object.isRequired,
  dimension: PropTypes.string.isRequired,
  onChange: PropTypes.func,
  disabled: PropTypes.bool
}
