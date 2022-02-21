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
import {
  TextField,
  makeStyles,
  Box,
  FormControlLabel, Checkbox, IconButton, InputAdornment, MenuItem
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import AutoComplete from '@material-ui/lab/Autocomplete'
import HelpOutlineIcon from '@material-ui/icons/HelpOutline'

const HelpButton = React.memo(function HelpButton(props) {
  return <IconButton size="small">
    {<HelpOutlineIcon fontSize='small'/>}
  </IconButton>
})

const useHelpAdornmentStyles = makeStyles(theme => ({
  root: {},
  withOtherAdornment: {
    marginRight: theme.spacing(3)
  }
}))

const HelpAdornment = React.memo(function HelpAdornment({withOtherAdornment}) {
  const classes = useHelpAdornmentStyles()
  return <InputAdornment
    position="end"
    className={withOtherAdornment ? classes.withOtherAdornment : classes.root}
  >
    <HelpButton />
  </InputAdornment>
})
HelpAdornment.propTypes = {
  withOtherAdornment: PropTypes.bool
}

const useWithHelpStyles = makeStyles(theme => ({
  root: {
    '&:not(:hover)': {
      '& #help': {
        display: 'none'
      }
    }
  }
}))

const TextFieldWithHelp = React.memo(function TextFieldWithHelp(props) {
  const {withOtherAdornment, ...otherProps} = props
  const classes = useWithHelpStyles()
  return <TextField
    className={classes.root}
    InputProps={{endAdornment: (
      <div id="help">
        <HelpAdornment withOtherAdornment={withOtherAdornment}/>
      </div>
    )}}
    {...otherProps}
  />
})
TextFieldWithHelp.propTypes = {
  withOtherAdornment: PropTypes.bool
}

const WithHelp = React.memo(function TextFieldWithHelp(props) {
  const classes = useWithHelpStyles()
  return <Box display="flex" alignItems="center" className={classes.root}>
    <Box flexGrow={1} {...props}/>
    <Box>
      <div id="help">
        <HelpButton />
      </div>
    </Box>
  </Box>
})

export const StringEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const label = otherProps.label || quantityDef.name
  const [value, setValue] = useState()

  useEffect(() => {
    setValue(section[quantityDef.name])
  }, [quantityDef, section])

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextFieldWithHelp fullWidth variant='filled' size='small'
    value={value || ''}
    placeholder={quantityDef.description}
    onChange={event => handleChange(event.target.value)} {...otherProps}
    label={label}
  />
})
StringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

const useNumberEditQuantityStyles = makeStyles(theme => ({
  unitSelect: {
    marginLeft: theme.spacing(1),
    width: '150px'
  }
}))

export const NumberEditQuantity = React.memo((props) => {
  const classes = useNumberEditQuantityStyles()
  const {quantityDef, section, onChange, minValue, maxValue, ...otherProps} = props
  const label = otherProps.label || quantityDef.name
  const [value, setValue] = useState()
  const [convertedValue, setConvertedValue] = useState()
  const [error, setError] = useState('')
  const systemUnits = useUnits()
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : '')
  const dimension = quantityDef.unit && unitMap[quantityDef.unit].dimension
  const units = quantityDef.unit && conversionMap[dimension].units
  const isUnit = quantityDef.unit && ['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)
  const timeout = useRef()

  useEffect(() => {
    let newValue = section[quantityDef.name]
    newValue = `${(isUnit ? (!isNaN(Number(newValue)) || newValue === '' ? convertUnit(Number(newValue), quantityDef.unit, unit) : '') : newValue)}`
    setValue(newValue)
    setConvertedValue((newValue === '' || newValue === undefined ? `${defaultValue}` : Number(newValue)))
  }, [defaultValue, isUnit, quantityDef, section, unit])

  const handleChangeUnit = useCallback((newUnit) => {
    setUnit(newUnit)
    setConvertedValue(`${(isUnit ? (!isNaN(Number(value)) || value === '' ? convertUnit(Number(value), quantityDef.unit, newUnit) : '') : value)}`)
  }, [isUnit, quantityDef, value])

  const isValidNumber = useCallback((value) => {
    if (['int64', 'int32', 'int'].includes(quantityDef.type?.type_data)) {
      const num = Number(value)
      return Number.isInteger(num)
    } else if (['uint64', 'uint32', 'uint'].includes(quantityDef.type?.type_data)) {
      const num = Number(value)
      return Number.isInteger(num) && num > 0
    } else if (['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)) {
      const num = Number(value)
      return !isNaN(num)
    }
  }, [quantityDef])

  const validation = useCallback((newValue) => {
    setError('')
    if (newValue === '') {
      if (defaultValue !== '') {
        setConvertedValue(`${(isUnit ? (!isNaN(Number(defaultValue)) || value === '' ? convertUnit(Number(defaultValue), unit, quantityDef.unit) : '') : defaultValue)}`)
        setValue(Number(defaultValue))
      } else {
        setConvertedValue(`${defaultValue}`)
        setValue(defaultValue)
      }
    } else if (!isValidNumber(newValue)) {
      setError('Please enter a valid number!')
    } else if (minValue !== undefined && Number(newValue) < minValue) {
      setError(`The value should be higher than or equal to ${minValue}`)
    } else if (maxValue !== undefined && Number(newValue) > maxValue) {
      setError(`The value should be less than or equal to ${maxValue}`)
    } else {
      setValue((isUnit ? (!isNaN(Number(newValue)) || value === '' ? convertUnit(Number(newValue), unit, quantityDef.unit) : '') : newValue))
      setConvertedValue(`${Number(newValue)}`)
    }
  }, [defaultValue, isUnit, isValidNumber, maxValue, minValue, quantityDef, unit, value])

  const handleChangeValue = useCallback((newValue) => {
    setConvertedValue(`${newValue}`)
    if (onChange) {
      onChange((newValue === '' ? Number(defaultValue) : (isUnit ? (!isNaN(Number(newValue)) || value === '' ? convertUnit(Number(newValue), unit, quantityDef.unit) : '') : newValue)), section, quantityDef)
    }
    clearTimeout(timeout.current)
    timeout.current = setTimeout(() => {
      validation(newValue)
    }, 1000)
  }, [isUnit, value, validation, unit, defaultValue, onChange, quantityDef, section, timeout])

  const handleValidator = useCallback((event) => {
    validation(event.target.value)
  }, [validation])

  return <Box display='flex'>
    <TextFieldWithHelp fullWidth variant='filled' size='small'
      value={convertedValue || ''}
      onBlur={handleValidator} error={!!error} helperText={error}
      placeholder={quantityDef.description}
      onChange={event => handleChangeValue(event.target.value)}
      {...otherProps} label={label}
    />
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
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const EnumEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const label = otherProps.label || quantityDef.name
  const [value, setValue] = useState(section[quantityDef.name] || quantityDef.default || '')

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value === '' ? undefined : value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextFieldWithHelp
    select variant='filled' size='small' withOtherAdornment fullWidth
    label={label} {...otherProps} value={value}
    onChange={event => handleChange(event.target.value)}
  >
    {quantityDef.type?.type_data.map(item => <MenuItem value={item} key={item}>{item}</MenuItem>)}
  </TextFieldWithHelp>
})
EnumEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const AutocompleteEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const label = otherProps.label || quantityDef.name
  const [value, setValue] = useState(section[quantityDef.name] || quantityDef.default || null)

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange((value === '' ? undefined : value), section, quantityDef)
    }
  }, [onChange, quantityDef, section, setValue])

  return <AutoComplete
    options={quantityDef.type.type_data}
    onChange={(event, value) => handleChange(value)}
    ListboxProps={{style: {maxHeight: '150px'}}}
    value={value}
    renderInput={params => (
      <TextField
        {...params}
        variant='filled' size='small' label={label}
        placeholder={quantityDef.description} fullWidth/>
    )}
    {...otherProps}
  />
})
AutocompleteEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const BoolEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const label = otherProps.label || quantityDef.name
  const [value, setValue] = useState()
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : '')

  useEffect(() => {
    setValue(section[quantityDef.name] || defaultValue)
  }, [defaultValue, quantityDef, section])

  const handleChange = useCallback((newValue) => {
    setValue(newValue)
    if (onChange) {
      onChange((newValue === '' ? defaultValue : newValue), section, quantityDef)
    }
  }, [defaultValue, onChange, quantityDef, section])

  return <WithHelp>
    <FormControlLabel
      label={label}
      control={<Checkbox onChange={event => handleChange(event.target.checked)} color="primary" checked={(!!value)}/>}
    />
  </WithHelp>
})
BoolEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
