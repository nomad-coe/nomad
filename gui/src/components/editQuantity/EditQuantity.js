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
import React, {useCallback, useEffect, useState} from 'react'
import {TextField, makeStyles, InputAdornment, Select, Box, FormControl, InputLabel} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit} from '../../units'
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
    marginLeft: 3,
    width: '100px'
  }
}))

export const StringEditQuantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, onChange} = props
  const [value, setValue] = useState()

  const eln = quantityDef?.m_annotations?.eln || []
  const label = eln[0]?.label || quantityDef.name
  const sectionProps = eln[0]?.props || undefined
  const multiline = sectionProps?.multiline
  const minRows = sectionProps?.minRows

  useEffect(() => {
    setValue(section[quantityDef.name])
  }, [quantityDef, section])

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextField fullWidth='true' variant='filled' size='small'
    multiline={multiline} minRows={minRows}
    value={value || ''}
    label={label}
    InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{quantityDef?.unit}</InputAdornment>}}
    placeholder={quantityDef?.description}
    onChange={event => handleChange(event.target.value)}>
  </TextField>
})
StringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const NumberEditQuantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, onChange} = props
  const [value, setValue] = useState()
  const [error, setError] = useState('')
  const [selectedUnit, setSelectedUnit] = useState(quantityDef?.unit)

  const eln = quantityDef?.m_annotations?.eln || []
  const sectionProps = eln[0]?.props || {}
  const label = eln[0]?.label || quantityDef.name
  const maxValue = sectionProps?.maxValue
  const minValue = sectionProps?.minValue
  const defaultValue = (sectionProps?.defaultValue !== undefined ? sectionProps?.defaultValue : '')
  const dimension = quantityDef?.unit && unitMap[quantityDef?.unit].dimension
  const units = quantityDef?.unit && conversionMap[dimension].units

  useEffect(() => {
    setValue(section[quantityDef.name] || `${defaultValue}`)
  }, [defaultValue, quantityDef, section])

  let timeout = null
  const handleChange = useCallback((value, unit) => {
    setValue(value)
    setSelectedUnit(unit)
    if (onChange) {
      onChange((value === '' ? defaultValue : (quantityDef?.unit ? convertUnit(Number(value), unit, quantityDef.unit) : value)), section, quantityDef)
    }
    clearTimeout(timeout)
    timeout = setTimeout(() => {
      validation(value)
    }, 1000)
  }, [defaultValue, onChange, quantityDef, section])

  const isValidNumber = useCallback((value) => {
    if (['int64', 'int32', 'int'].includes(quantityDef?.type?.type_data)) {
      const num = Number(value)
      return Number.isInteger(num)
    } else if (['uint64', 'uint32', 'uint'].includes(quantityDef?.type?.type_data)) {
      const num = Number(value)
      return Number.isInteger(num) && num > 0
    } else if (['float64', 'float32', 'float'].includes(quantityDef?.type?.type_data)) {
      const num = Number(value)
      return !isNaN(num)
    }
  }, [quantityDef])

  const validation = useCallback((value) => {
    setError('')
    if (value === '') {
      setValue(`${defaultValue}`)
    } else if (!isValidNumber(value)) {
      setError('Please enter a valid number!')
    } else if (minValue !== undefined && Number(value) < minValue) {
      setError(`The value should be higher than or equal to ${minValue}`)
    } else if (maxValue !== undefined && Number(value) > maxValue) {
      setError(`The value should be less than or equal to ${maxValue}`)
    }
  }, [defaultValue, isValidNumber, maxValue, minValue])

  const handleValidator = useCallback((event) => {
    validation(event.target.value)
  }, [validation])

  return <Box display='flex'>
    <TextField fullWidth='true' variant='filled' size='small'
      value={value || ''}
      label={label}
      onBlur={handleValidator} error={!!error} helperText={error}
      placeholder={quantityDef?.description}
      onChange={event => handleChange(event.target.value, selectedUnit)}>
    </TextField>
    {quantityDef?.unit && <FormControl className={classes.unitSelect} variant='filled' size='small'>
      <InputLabel htmlFor={'unit'} shrink={true}>{'unit'}</InputLabel>
      <Select native value={selectedUnit}
        onChange={(event) => handleChange(value, event.target.value)}>
        {units.map(unit => <option key={unit}>{(new Unit(unit)).label()}</option>)}
      </Select>
    </FormControl>}
  </Box>
})
NumberEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const EnumEditQuantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, onChange} = props
  const [value, setValue] = useState()
  const [selectedUnit, setSelectedUnit] = useState(quantityDef?.unit)

  const eln = quantityDef?.m_annotations?.eln || []
  const sectionProps = eln[0]?.props || {}
  const label = eln[0]?.label || quantityDef.name
  const defaultValue = (sectionProps?.defaultValue !== undefined ? sectionProps?.defaultValue : '')
  const dimension = quantityDef?.unit && unitMap[quantityDef?.unit].dimension
  const units = quantityDef?.unit && conversionMap[dimension].units

  useEffect(() => {
    setValue(section[quantityDef.name] || defaultValue)
  }, [defaultValue, quantityDef, section])

  const handleChange = useCallback((value, unit) => {
    setValue(value)
    setSelectedUnit(unit)
    if (onChange) {
      onChange((value === '' ? defaultValue : (quantityDef?.unit ? convertUnit(Number(value), unit, quantityDef.unit) : value)), section, quantityDef)
    }
  }, [defaultValue, onChange, quantityDef, section])

  return <Box display='flex'>
    <FormControl variant='filled' size='small' fullWidth>
      <InputLabel htmlFor={label} shrink={value !== undefined && value !== ''}>{label}</InputLabel>
      <Select native value={value}
        onChange={event => handleChange(event.target.value, selectedUnit)}>
        {defaultValue === '' && <option key={''}>{''}</option>}
        {quantityDef?.type?.type_data.map(item => <option key={item}>{item}</option>)}
      </Select>
    </FormControl>
    {quantityDef?.unit && <FormControl className={classes.unitSelect} variant='filled' size='small'>
      <InputLabel htmlFor={'unit'} shrink={true}>{'unit'}</InputLabel>
      <Select native value={selectedUnit}
        onChange={(event) => handleChange(value, event.target.value)}>
        {units.map(unit => <option key={unit}>{(new Unit(unit)).label()}</option>)}
      </Select>
    </FormControl>}
  </Box>
})
EnumEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const EditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange} = props

  const eln = quantityDef?.m_annotations?.eln
  const component = (eln.length > 0 ? eln[0]?.component : undefined)

  if (component === 'StringEditQuantity') {
    return <StringEditQuantity quantityDef={quantityDef} section={section} onChange={onChange}/>
  } else if (component === 'NumberEditQuantity') {
    return <NumberEditQuantity quantityDef={quantityDef} section={section} onChange={onChange}/>
  } else if (component === 'EnumEditQuantity') {
    return <EnumEditQuantity quantityDef={quantityDef} section={section} onChange={onChange}/>
  }
})
EditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
