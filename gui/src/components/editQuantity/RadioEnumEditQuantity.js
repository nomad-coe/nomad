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
  FormControl,
  FormLabel, RadioGroup, Radio
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {getFieldProps, WithHelp} from './StringEditQuantity'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const RadioEnumEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, ...otherProps} = props
  const fieldProps = getFieldProps(quantityDef)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value === '' ? undefined : value)
    }
  }, [onChange])

  return (
    <WithHelp {...fieldProps} label={label}>
      <FormControl>
        <FormLabel>{label}</FormLabel>
        <RadioGroup row>
          {quantityDef.type?.type_data.map(item => (
            <FormControlLabel
              value={item} key={item}
              control={(
                <Radio
                  checked={value === item}
                  onClick={() => handleChange(item)}
                  {...otherProps}
                />
              )}
              label={item}
            />))}
        </RadioGroup>
      </FormControl>
    </WithHelp>
  )
})
RadioEnumEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}
