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
import React, {useCallback, useEffect, useMemo, useState} from 'react'
import {TextField, makeStyles, Box, MenuItem} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import {debounce} from 'lodash'
import {TextFieldWithHelp, getFieldProps} from './StringEditQuantity'

export const useNumberEditQuantityStyles = makeStyles(theme => ({
  unitSelect: {
    marginLeft: theme.spacing(1),
    width: '150px'
  }
}))

export const NumberFieldWithUnit = React.memo((props) => {
  const {onChange, defaultUnit, dataType, minValue, maxValue, unit, defaultValue, ...otherProps} = props
  const [convertedValue, setConvertedValue] = useState()
  const [error, setError] = useState('')
  const isUnit = unit !== undefined

  useEffect(() => {
    if (defaultValue === undefined || defaultValue === '' || isNaN(Number(defaultValue))) {
      setConvertedValue('')
    } else {
      setConvertedValue((isUnit ? convertUnit(Number(defaultValue), defaultUnit, unit) : defaultValue))
    }
  }, [defaultValue, isUnit, defaultUnit, unit])

  const isValidNumber = useCallback((value) => {
    if (['int64', 'int32', 'int'].includes(dataType)) {
      const num = Number(value)
      return Number.isInteger(num)
    } else if (['uint64', 'uint32', 'uint'].includes(dataType)) {
      const num = Number(value)
      return Number.isInteger(num) && num >= 0
    } else if (['float64', 'float32', 'float'].includes(dataType)) {
      const num = Number(value)
      return !isNaN(num)
    }
  }, [dataType])

  const validation = useCallback((val, fastEvaluation) => {
    setError('')
    let newValue = val.replace(/,/g, '.')
    if (newValue === '') {
      setConvertedValue('')
      if (onChange) onChange('')
    } else if (fastEvaluation) {
      if (!newValue.match(/^[+-]?((\d+|\.\d?|\d+\.|\d+\.\d+)|(\d+|\.\d?|\d+\.|\d+\.\d+)(e|e\+|e-)\d*)?$/)) setError('Please enter a valid number!')
    } else if (!isValidNumber(newValue)) {
      setError('Please enter a valid number!')
    } else {
      let originalValue = (isUnit ? convertUnit(Number(newValue), unit, defaultUnit) : newValue)
      if (minValue !== undefined && originalValue < minValue) {
        setError(`The value should be higher than or equal to ${minValue}${(isUnit ? `${(new Unit(defaultUnit)).label()}` : '')}`)
      } else if (maxValue !== undefined && originalValue > maxValue) {
        setError(`The value should be less than or equal to ${maxValue}${(isUnit ? `${(new Unit(defaultUnit)).label()}` : '')}`)
      } else {
        setConvertedValue(Number(newValue))
        if (onChange) onChange(originalValue)
      }
    }
  }, [isUnit, isValidNumber, maxValue, minValue, onChange, defaultUnit, unit])

  const debouncedValidation = useMemo(() => {
    return debounce(validation, 2000)
  }, [validation])

  const handleChangeValue = useCallback((newValue) => {
    setConvertedValue(newValue)
    validation(newValue, true)
    debouncedValidation(newValue)
  }, [debouncedValidation, validation])

  const handleValidator = useCallback((event) => {
    debouncedValidation.flush()
  }, [debouncedValidation])

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    value={convertedValue !== undefined ? convertedValue : ''}
    onBlur={handleValidator} error={!!error} helperText={error}
    onChange={event => handleChangeValue(event.target.value)}
    {...otherProps}
  />
})
NumberFieldWithUnit.propTypes = {
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  onChange: PropTypes.func.isRequired,
  defaultUnit: PropTypes.string,
  dataType: PropTypes.string,
  defaultValue: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
  unit: PropTypes.string
}

export const NumberEditQuantity = React.memo((props) => {
  const classes = useNumberEditQuantityStyles()
  const {quantityDef, section, onChange, ...otherProps} = props
  const systemUnits = useUnits()
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : '')
  const dimension = quantityDef.unit && unitMap[quantityDef.unit].dimension
  const units = quantityDef.unit && conversionMap[dimension].units
  const isUnit = quantityDef.unit && ['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)

  const handleChangeUnit = useCallback((newUnit) => {
    setUnit(newUnit)
  }, [])

  const handleChangeValue = useCallback((newValue) => {
    if (onChange) {
      onChange(newValue, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <Box display='flex'>
    <NumberFieldWithUnit
      onChange={handleChangeValue}
      defaultUnit={quantityDef.unit}
      dataType={quantityDef.type?.type_data}
      unit={unit}
      defaultValue={section[quantityDef.name] !== undefined ? section[quantityDef.name] : defaultValue}
      helpDescription={quantityDef.description}
      {...getFieldProps(quantityDef)}
      {...otherProps}/>
    {isUnit && <TextField
      className={classes.unitSelect} variant='filled' size='small' select
      label="unit" value={unit}
      onChange={(event) => handleChangeUnit(event.target.value)}
    >
      {units.map(unit => <MenuItem key={unit} value={unit}>{(new Unit(unit)).label()}</MenuItem>)}
    </TextField>}
  </Box>
})
NumberEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
