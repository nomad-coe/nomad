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
  FormControlLabel,
  Checkbox,
  IconButton,
  InputAdornment,
  MenuItem,
  Dialog,
  DialogContent,
  FormControl,
  FormLabel, RadioGroup, Radio, Slider, DialogTitle, Collapse
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import AutoComplete from '@material-ui/lab/Autocomplete'
import HelpOutlineIcon from '@material-ui/icons/HelpOutline'
import DialogActions from '@material-ui/core/DialogActions'
import Button from '@material-ui/core/Button'
import Markdown from '../Markdown'
import {dateFormat} from '../../config'
import {KeyboardDatePicker, KeyboardTimePicker} from '@material-ui/pickers'
import {getTime} from 'date-fns'
import AccessTimeIcon from '@material-ui/icons/AccessTime'
import ArrowDownIcon from '@material-ui/icons/ArrowDropDown'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'

const HelpDialog = React.memo(({title, description}) => {
  const [open, setOpen] = useState(false)

  return <React.Fragment>
    {description && <IconButton size="small" onClick={() => setOpen(true)}>
      {<HelpOutlineIcon fontSize='small'/>}
    </IconButton>}
    {open && <Dialog open={open}>
      <DialogTitle>
        {title}
      </DialogTitle>
      <DialogContent>
        <Markdown>{description}</Markdown>
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={() => setOpen(false)} color="secondary">
          Close
        </Button>
      </DialogActions>

    </Dialog>}
  </React.Fragment>
})
HelpDialog.propTypes = {
  title: PropTypes.string,
  description: PropTypes.string
}

const useHelpAdornmentStyles = makeStyles(theme => ({
  root: {},
  withOtherAdornment: {
    marginRight: theme.spacing(3)
  }
}))

const HelpAdornment = React.memo(({title, description, withOtherAdornment}) => {
  const classes = useHelpAdornmentStyles()
  return <InputAdornment
    position="end"
    className={withOtherAdornment ? classes.withOtherAdornment : classes.root}
  >
    <HelpDialog title={title} description={description}/>
  </InputAdornment>
})
HelpAdornment.propTypes = {
  withOtherAdornment: PropTypes.bool,
  title: PropTypes.string,
  description: PropTypes.string
}

const useWithHelpStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    '&:not(:hover)': {
      '& #help': {
        display: 'none'
      }
    }
  }
}))

const TextFieldWithHelp = React.memo((props) => {
  const {withOtherAdornment, label, helpDescription, ...otherProps} = props
  const classes = useWithHelpStyles()
  return <TextField
    className={classes.root}
    InputProps={(helpDescription && {endAdornment: (
      <div id="help">
        <HelpAdornment title={label} description={helpDescription} withOtherAdornment={withOtherAdornment}/>
      </div>
    )})}
    label={label}
    {...otherProps}
  />
})
TextFieldWithHelp.propTypes = {
  withOtherAdornment: PropTypes.bool,
  label: PropTypes.string,
  helpDescription: PropTypes.string
}

const WithHelp = React.memo((props) => {
  const {label, helpDescription, ...otherProps} = props
  const classes = useWithHelpStyles()
  if (!helpDescription) {
    return ''
  }
  return <Box display="flex" alignItems="center" className={classes.root}>
    <Box flexGrow={1} {...otherProps}/>
    <Box>
      <div id="help">
        <HelpDialog title={label} description={helpDescription} />
      </div>
    </Box>
  </Box>
})
WithHelp.propTypes = {
  label: PropTypes.string,
  helpDescription: PropTypes.string
}

const capitalize = (s) => {
  if (typeof s !== 'string') return ''
  return s.charAt(0).toUpperCase() + s.slice(1)
}

function getArchiveValue(quantityDef, section) {
  let value = section[quantityDef.name]
  if (value === undefined) {
    return quantityDef.default
  }
  return value
}

function getFieldProps(quantityDef) {
  const eln = quantityDef?.m_annotations?.eln
  let name = quantityDef.name.replace(/_/g, ' ')
  let capitalizeName = capitalize(name)
  let label = (eln.length > 0 ? eln[0]?.label : undefined) || capitalizeName
  return {
    label: label,
    helpDescription: quantityDef.description
  }
}

export const StringEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const value = getArchiveValue(quantityDef, section)

  const handleChange = useCallback((newValue) => {
    if (onChange) {
      onChange(newValue === '' ? undefined : newValue, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    defaultValue={value !== undefined ? value : ''}
    onChange={event => handleChange(event.target.value)}
    {...getFieldProps(quantityDef)}
    {...otherProps}
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
      minValue={minValue}
      maxValue={maxValue}
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
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

const NumberFieldWithUnit = React.memo((props) => {
  const {onChange, defaultUnit, dataType, minValue, maxValue, unit, defaultValue, ...otherProps} = props
  const [convertedValue, setConvertedValue] = useState()
  const [error, setError] = useState('')
  const timeout = useRef()
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
      if (!newValue.match(/^[+-]?((\d+|\.\d?|\d+\.|\d+\.\d+)|(\d+|\.\d?|\d+\.|\d+\.\d+)(e|e\+|e-)\d?)?$/)) setError('Please enter a valid number!')
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

  const handleChangeValue = useCallback((newValue) => {
    setConvertedValue(newValue)
    validation(newValue, true)
    clearTimeout(timeout.current)
    timeout.current = setTimeout(() => {
      validation(newValue, false)
    }, 3000)
  }, [validation])

  const handleValidator = useCallback((event) => {
    clearTimeout(timeout.current)
    validation(event.target.value, false)
  }, [validation])

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

const StringField = React.memo((props) => {
  const {onChange, defaultValue, ...otherProps} = props

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    onChange={event => onChange(event.target.value)}
    defaultValue={defaultValue}
    {...otherProps}
  />
})
StringField.propTypes = {
  onChange: PropTypes.func.isRequired,
  defaultValue: PropTypes.string
}

export const EnumEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const [value, setValue] = useState(section[quantityDef.name] || quantityDef.default || '')

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value === '' ? undefined : value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextFieldWithHelp
    select variant='filled' size='small' withOtherAdornment fullWidth
    value={value}
    onChange={event => handleChange(event.target.value)}
    {...getFieldProps(quantityDef)}
    {...otherProps}
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
      <TextFieldWithHelp
        {...params}
        variant='filled' size='small'
        placeholder={quantityDef.description} fullWidth
        {...getFieldProps(quantityDef)}
        {...otherProps}
      />
    )}
  />
})
AutocompleteEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const RadioButtonEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const [value, setValue] = useState(section[quantityDef.name] || quantityDef.default || '')

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value === '' ? undefined : value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <FormControl>
    <FormLabel>{getFieldProps(quantityDef).label}</FormLabel>
    <RadioGroup row>
      {quantityDef.type?.type_data.map(item => <FormControlLabel value={item} key={item} control={<Radio checked={value === item} onClick={event => handleChange(item)} {...otherProps}/>} label={item}/>)}
    </RadioGroup>
  </FormControl>
})
RadioButtonEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const BoolEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
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

  return <WithHelp {...getFieldProps(quantityDef)}>
    <FormControlLabel
      control={<Checkbox onChange={event => handleChange(event.target.checked)} color="primary" checked={(!!value)} {...otherProps}/>}
      label={getFieldProps(quantityDef).label}
    />
  </WithHelp>
})
BoolEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const SliderEditQuantity = React.memo((props) => {
  const classes = useNumberEditQuantityStyles()
  const {quantityDef, section, onChange, minValue, maxValue, ...otherProps} = props
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : undefined)
  const [value, setValue] = useState(0)
  const [convertedValue, setConvertedValue] = useState(0)
  const dimension = quantityDef.unit && unitMap[quantityDef.unit].dimension
  const units = quantityDef.unit && conversionMap[dimension].units
  const systemUnits = useUnits()
  const isUnit = quantityDef.unit && ['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)

  useEffect(() => {
    let newValue = section[quantityDef.name] !== undefined ? section[quantityDef.name] : (defaultValue || minValue)
    setValue(newValue)
    setConvertedValue(`${(isUnit ? (!isNaN(Number(newValue)) || newValue === '' ? convertUnit(Number(newValue), quantityDef.unit, unit) : '') : newValue)}`)
  }, [defaultValue, isUnit, minValue, quantityDef, section, unit])

  const handleChangeUnit = useCallback((newUnit) => {
    setUnit(newUnit)
    setConvertedValue(`${(isUnit ? (!isNaN(Number(value)) || value === '' ? convertUnit(Number(value), quantityDef.unit, newUnit) : '') : value)}`)
  }, [isUnit, quantityDef, value])

  const handleChangeValue = useCallback((event, newValue) => {
    if (typeof newValue !== 'number') return
    setConvertedValue(`${newValue}`)
    if (onChange) {
      onChange((isUnit ? (newValue === '' ? newValue : (!isNaN(Number(newValue)) ? convertUnit(Number(newValue), unit, quantityDef.unit) : '')) : newValue), section, quantityDef)
    }
    setValue((isUnit ? (!isNaN(Number(newValue)) || newValue === '' ? convertUnit(Number(newValue), unit, quantityDef.unit) : '') : newValue))
    setConvertedValue(`${Number(newValue)}`)
  }, [isUnit, unit, onChange, quantityDef, section])

  return <FormControl fullWidth>
    <FormLabel>{getFieldProps(quantityDef).label}</FormLabel>
    <Box display='flex'>
      <Slider
        value={Number(convertedValue)}
        min={convertUnit(Number(minValue), quantityDef.unit, unit)}
        max={convertUnit(Number(maxValue), quantityDef.unit, unit)}
        onChange={handleChangeValue}
        valueLabelDisplay={(!isUnit ? 'on' : 'off')}
        {...otherProps}/>
      {isUnit && <TextField
        className={classes.unitSelect} variant='filled' size='small' select
        label="unit" value={unit}
        onChange={(event) => handleChangeUnit(event.target.value)}
      >
        {units.map(unit => <MenuItem key={unit} value={unit}>{(new Unit(unit)).label()}</MenuItem>)}
      </TextField>}
    </Box>
  </FormControl>
})
SliderEditQuantity.propTypes = {
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

const useDatesEditQuantityStyles = makeStyles(theme => ({
  startDate: {
  },
  endDate: {
    marginLeft: theme.spacing(1)
  }
}))

export const DateTimeEditQuantity = React.memo((props) => {
  const classes = useDatesEditQuantityStyles()
  const {quantityDef, section, onChange, format, time, ...otherProps} = props
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : '')
  const [value, setValue] = useState()
  const [current, setCurrent] = useState()
  const [error, setError] = useState('')

  useEffect(() => {
    setValue(section[quantityDef.name] || defaultValue || null)
  }, [defaultValue, quantityDef, section])

  const handleAccept = useCallback((newValue) => {
    if (newValue !== null && newValue !== undefined && isNaN(getTime(newValue))) {
      setError('Invalid date format.')
      return
    }
    setError('')
    if (newValue !== undefined) setValue(newValue)
    if (onChange) {
      onChange(newValue || '', section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  const handleChange = useCallback((newValue) => {
    setCurrent(newValue)
  }, [])

  const handleBlur = useCallback(() => {
    handleAccept(current)
  }, [current, handleAccept])

  const renderProps = {
    className: classes.startDate,
    size: 'small',
    error: !!error,
    variant: 'inline',
    inputVariant: 'filled',
    fullWidth: true,
    label: getFieldProps(quantityDef).label,
    value: value,
    invalidDateMessage: error,
    onAccept: handleAccept,
    onChange: handleChange,
    onBlur: handleBlur,
    onKeyDown: (event) => { if (event.key === 'Enter') { handleAccept(current) } },
    ...otherProps
  }

  if (time) {
    return <KeyboardTimePicker
      {...renderProps}
      format={format || `HH:mm`}
      keyboardIcon={<AccessTimeIcon />}
    />
  } else {
    return <KeyboardDatePicker
      {...renderProps}
      format={format || `${dateFormat} HH:mm`}
    />
  }
})
DateTimeEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired,
  format: PropTypes.string,
  time: PropTypes.bool
}

export const DateEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props

  return <DateTimeEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} format={dateFormat} {...otherProps}/>
})
DateEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const TimeEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props

  return <DateTimeEditQuantity quantityDef={quantityDef} section={section} onChange={onChange} time {...otherProps}/>
})
TimeEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

const ListEditQuantity = React.memo((props) => {
  const {quantityDef, section, component, componentProps, defaultValues, collapse, onChange, actions} = props
  const [values, setValues] = useState(defaultValues)
  const [open, setOpen] = useState((collapse !== undefined ? !collapse : false))
  let Component = component

  const handleChange = useCallback((newValue, index) => {
    let newValues = [...values]
    newValues[index] = newValue
    setValues(newValues)
    if (onChange) {
      onChange(newValues, section, quantityDef)
    }
  }, [onChange, quantityDef, section, values])

  if (!defaultValues) return ''

  return <Box display='block'>
    <Box display='flex'>
      <WithHelp helpDescription={quantityDef.description} {...getFieldProps(quantityDef)}>
        <Box sx={{ flexDirection: 'column' }} onClick={() => setOpen(!open)} height={'26px'}>
          <FormControlLabel
            control={open ? <ArrowDownIcon/> : <ArrowRightIcon/>}
            label={getFieldProps(quantityDef).label}
          />
        </Box>
      </WithHelp>
      {open && actions}
    </Box>
    <Collapse in={open}>
      <div>
        {defaultValues.map((value, index) => {
          return <Box key={index} marginTop={1}>
            <Component
              defaultValue={value !== undefined ? value : ''} onChange={newValue => handleChange(newValue, index)} {...componentProps} label={undefined}
              inputProps={{style: { padding: '14px' }}}
            />
          </Box>
        })}
      </div>
    </Collapse>
  </Box>
})
ListEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  component: PropTypes.any.isRequired,
  componentProps: PropTypes.object.isRequired,
  defaultValues: PropTypes.arrayOf(PropTypes.any),
  onChange: PropTypes.func.isRequired,
  collapse: PropTypes.bool,
  actions: PropTypes.element
}

export const ListNumberEditQuantity = React.memo((props) => {
  const classes = useNumberEditQuantityStyles()
  const {quantityDef, section, onChange, minValue, maxValue, ...otherProps} = props
  const systemUnits = useUnits()
  const shape = quantityDef.type?.shape
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : Array.apply(null, Array(shape[0])).map(() => ''))
  const dimension = quantityDef.unit && unitMap[quantityDef.unit].dimension
  const units = quantityDef.unit && conversionMap[dimension].units
  const isUnit = quantityDef.unit && ['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)
  let values = section[quantityDef.name] || defaultValue

  const handleChangeUnit = useCallback((newUnit) => {
    setUnit(newUnit)
  }, [])

  const componentProps = {
    defaultUnit: quantityDef.unit,
    dataType: quantityDef.type?.type_data,
    minValue: minValue,
    maxValue: maxValue,
    unit: unit,
    ...otherProps
  }

  return <ListEditQuantity
    quantityDef={quantityDef}
    section={section}
    component={NumberFieldWithUnit}
    componentProps={componentProps}
    defaultValues={values}
    onChange={onChange}
    actions={isUnit && <TextField
      className={classes.unitSelect} variant='filled' size='small' select
      label="unit" value={unit}
      onChange={(event) => handleChangeUnit(event.target.value)}
    >
      {units.map(unit => <MenuItem key={unit} value={unit}>{(new Unit(unit)).label()}</MenuItem>)}
    </TextField>}/>
})
ListNumberEditQuantity.propTypes = {
  maxValue: PropTypes.number,
  minValue: PropTypes.number,
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const ListStringEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const shape = quantityDef.type?.shape
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : Array.apply(null, Array(shape[0])).map(() => ''))
  let values = section[quantityDef.name] || defaultValue

  return <ListEditQuantity
    quantityDef={quantityDef}
    section={section}
    component={StringField}
    componentProps={otherProps}
    defaultValues={values}
    onChange={onChange}/>
})
ListStringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
