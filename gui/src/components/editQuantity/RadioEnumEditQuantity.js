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
import React, {useCallback, useState} from 'react'
import {
  FormControlLabel,
  FormControl,
  FormLabel, RadioGroup, Radio, Box
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {getFieldProps} from './StringEditQuantity'

export const RadioEnumEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const [value, setValue] = useState(section[quantityDef.name] || quantityDef.default || '')

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value === '' ? undefined : value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return (
    <Box marginTop={2}>
      <FormControl>
        <FormLabel>{getFieldProps(quantityDef).label}</FormLabel>
        <RadioGroup row>
          {quantityDef.type?.type_data.map(item => <FormControlLabel value={item} key={item} control={<Radio checked={value === item} onClick={event => handleChange(item)} {...otherProps}/>} label={item}/>)}
        </RadioGroup>
      </FormControl>
    </Box>
  )
})
RadioEnumEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
