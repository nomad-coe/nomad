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
import React, {useCallback, useEffect, useRef, useState} from 'react'
import {TextField, makeStyles, InputAdornment, Select, Box, FormControl} from '@material-ui/core'
import PropTypes from 'prop-types'
import {toUnitSystem, Unit} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'

const useStyles = makeStyles(theme => ({
  editQuantity: {
    display: 'block',
    width: '100%'
  },
  adornment: {
    marginRight: theme.spacing(3)
  },
  unitSelect: {
    marginLeft: 3
  }
}))

export const StringEditQuantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, onChange, multiline, minRows} = props
  const [value, setValue] = useState()

  useEffect(() => {
    setValue(section[quantityDef.key])
  }, [quantityDef.key, section])

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextField fullWidth='true' variant='filled' size='small'
    multiline={multiline} minRows={minRows}
    value={value || ''}
    label={quantityDef?.name}
    InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{quantityDef?.unit}</InputAdornment>}}
    placeholder={quantityDef?.description}
    onChange={event => handleChange(event.target.value)}>
  </TextField>
})
StringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired,
  multiline: PropTypes.bool,
  minRows: PropTypes.number
}
StringEditQuantity.defaultProps = {
  multiline: false,
  minRows: 4
}

export const FloatEditQuantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, onChange, defaultValue, minValue, maxValue} = props
  const [value, setValue] = useState()
  const [error, setError] = useState('')
  const [selectedUnit, setSelectedUnit] = useState(quantityDef?.unit)

  const dimension = quantityDef?.unit && unitMap[quantityDef?.unit].dimension
  const units = quantityDef?.unit && conversionMap[dimension].units

  useEffect(() => {
    setValue(section[quantityDef.key])
  }, [quantityDef.key, section])

  const isFloat = useCallback((value) => {
    const num = Number(value)
    return !isNaN(num)
  }, [])

  const validation = useCallback((value) => {
    setError('')
    if (value === '') {
      setValue(`${defaultValue}`)
    } else if (!isFloat(value)) {
      setError('Please enter a valid number!')
    } else if (minValue !== undefined && Number(value) < minValue) {
      setError(`The value should be higher than ${minValue}`)
    } else if (maxValue !== undefined && Number(value) > maxValue) {
      setError(`The value should be less than ${maxValue}`)
    }
  }, [defaultValue, isFloat, maxValue, minValue])

  const timeout = useRef(null)
  const handleChange = useCallback((value, unit) => {
    setValue(value)
    setSelectedUnit(unit)
    if (onChange) {
      onChange((value === '' ? defaultValue : (quantityDef?.unit ? toUnitSystem(Number(value), unit, quantityDef.unit) : value)), section, quantityDef)
    }
    if (timeout.current) {
      clearTimeout(timeout.current)
    }
    timeout.current = setTimeout(() => {
      validation(value)
    }, 1000)
  }, [onChange, quantityDef, section, timeout, defaultValue, validation])

  const handleFloatValidator = useCallback((event) => {
    validation(event.target.value)
  }, [validation])

  return <Box display='flex'>
    <TextField fullWidth='true' variant='filled' size='small'
      value={value || ''}
      label={quantityDef?.name}
      onBlur={handleFloatValidator} error={!!error} helperText={error}
      placeholder={quantityDef?.description}
      onChange={event => handleChange(event.target.value, selectedUnit)}>
    </TextField>
    {quantityDef?.unit && <FormControl className={classes.unitSelect} variant='filled' size='small'>
      <Select native value={selectedUnit}
        onChange={(event) => handleChange(value, event.target.value)}>
        {units.map(unit => <option key={unit}>{(new Unit(unit)).label()}</option>)}
      </Select>
    </FormControl>}
  </Box>
})
FloatEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired,
  defaultValue: PropTypes.object,
  minValue: PropTypes.number,
  maxValue: PropTypes.number
}
FloatEditQuantity.defaultProps = {
  defaultValue: ''
}

export const IntegerEditQuantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, onChange, defaultValue, minValue, maxValue} = props
  const [value, setValue] = useState()
  const [error, setError] = useState('')
  const [selectedUnit, setSelectedUnit] = useState(quantityDef?.unit)

  const dimension = quantityDef?.unit && unitMap[quantityDef?.unit].dimension
  const units = quantityDef?.unit && conversionMap[dimension].units

  useEffect(() => {
    setValue(section[quantityDef.key])
  }, [quantityDef.key, section])

  const isInteger = useCallback((value) => {
    const num = Number(value)
    return Number.isInteger(num)
  }, [])

  const validation = useCallback((value) => {
    setError('')
    if (value === '') {
      setValue(`${defaultValue}`)
    } else if (!isInteger(value)) {
      setError('Please enter a valid number!')
    } else if (minValue !== undefined && Number(value) < minValue) {
      setError(`The value should be higher than ${minValue}`)
    } else if (maxValue !== undefined && Number(value) > maxValue) {
      setError(`The value should be less than ${maxValue}`)
    }
  }, [defaultValue, isInteger, maxValue, minValue])

  const timeout = useRef(null)
  const handleChange = useCallback((value, unit) => {
    setValue(value)
    setSelectedUnit(unit)
    if (onChange) {
      onChange((value === '' ? defaultValue : (quantityDef?.unit ? toUnitSystem(Number(value), unit, quantityDef.unit) : value)), section, quantityDef)
    }
    if (timeout.current) {
      clearTimeout(timeout.current)
    }
    timeout.current = setTimeout(() => {
      validation(value)
    }, 1000)
  }, [onChange, quantityDef, section, timeout, defaultValue, validation])

  const handleFloatValidator = useCallback((event) => {
    validation(event.target.value)
  }, [validation])

  return <Box display='flex'>
    <TextField fullWidth='true' variant='filled' size='small'
      value={value || ''}
      label={quantityDef?.name}
      onBlur={handleFloatValidator} error={!!error} helperText={error}
      placeholder={quantityDef?.description}
      onChange={event => handleChange(event.target.value, selectedUnit)}>
    </TextField>
    {quantityDef?.unit && <FormControl className={classes.unitSelect} variant='filled' size='small'>
      <Select native value={selectedUnit}
        onChange={(event) => handleChange(value, event.target.value)}>
        {units.map(unit => <option key={unit}>{(new Unit(unit)).label()}</option>)}
      </Select>
    </FormControl>}
  </Box>
})
IntegerEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired,
  defaultValue: PropTypes.object,
  minValue: PropTypes.number,
  maxValue: PropTypes.number
}
IntegerEditQuantity.defaultProps = {
  defaultValue: ''
}
