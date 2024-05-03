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
import PropTypes from 'prop-types'
import {dateFormat} from '../../config'
import {KeyboardDatePicker, KeyboardTimePicker} from '@material-ui/pickers'
import AccessTimeIcon from '@material-ui/icons/AccessTime'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const DateTimeEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, format, time, ...otherProps} = props
  const [dateValue, setDateValue] = useState(value || null)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  useEffect(() => {
    setDateValue(date => value
      ? (date?.toJSON?.() || date) === value ? date : value
      : null
    )
  }, [value])

  const handleChange = useCallback((date) => {
    setDateValue(date)
    if (onChange) {
      if (date === null) {
        onChange(undefined)
      } else {
        const jsonValue = date?.toJSON()
        if (jsonValue) {
          onChange(jsonValue)
        }
      }
    }
  }, [onChange, setDateValue])

  const renderProps = {
    size: 'small',
    variant: 'inline',
    inputVariant: 'filled',
    fullWidth: true,
    label: label,
    value: dateValue,
    onChange: handleChange,
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
  value: PropTypes.string,
  onChange: PropTypes.func,
  format: PropTypes.string,
  time: PropTypes.bool
}

export const DateEditQuantity = React.memo((props) => {
  return <DateTimeEditQuantity {...props} format={dateFormat}/>
})
DateEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}

export const TimeEditQuantity = React.memo((props) => {
  return <DateTimeEditQuantity time {...props}/>
})
TimeEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}
