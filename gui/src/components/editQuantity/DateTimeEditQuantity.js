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
import {makeStyles} from '@material-ui/core'
import PropTypes from 'prop-types'
import {dateFormat} from '../../config'
import {KeyboardDatePicker, KeyboardTimePicker} from '@material-ui/pickers'
import {getTime} from 'date-fns'
import AccessTimeIcon from '@material-ui/icons/AccessTime'
import {getFieldProps} from './StringEditQuantity'

const useDateTimeEditQuantityStyles = makeStyles(theme => ({
  startDate: {
  },
  endDate: {
    marginLeft: theme.spacing(1)
  }
}))

export const DateTimeEditQuantity = React.memo((props) => {
  const classes = useDateTimeEditQuantityStyles()
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
