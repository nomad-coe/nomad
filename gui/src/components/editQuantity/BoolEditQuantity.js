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
import React, {useCallback} from 'react'
import {
  FormControlLabel,
  Checkbox
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {getFieldProps, WithHelp} from './StringEditQuantity'

export const BoolEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : '')

  const handleChange = useCallback((newValue) => {
    if (onChange) {
      onChange((newValue === '' ? defaultValue : newValue), section, quantityDef)
    }
  }, [defaultValue, onChange, quantityDef, section])

  return <WithHelp {...getFieldProps(quantityDef)}>
    <FormControlLabel
      control={<Checkbox onChange={event => handleChange(event.target.checked)} color="primary" defaultChecked={!!section[quantityDef.name] || !!defaultValue} {...otherProps}/>}
      label={getFieldProps(quantityDef).label}
    />
  </WithHelp>
})
BoolEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
