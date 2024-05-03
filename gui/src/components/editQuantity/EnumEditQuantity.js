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
import {MenuItem} from '@material-ui/core'
import PropTypes from 'prop-types'
import {getFieldProps, TextFieldWithHelp} from './StringEditQuantity'
import AutoComplete from '@material-ui/lab/Autocomplete'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const EnumEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, suggestions, ...otherProps} = props
  const fieldProps = getFieldProps(quantityDef)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleChange = useCallback(value => {
    if (onChange) {
      onChange(value === '' ? undefined : value)
    }
  }, [onChange])

  if ((quantityDef.type?.type_kind === 'python' && quantityDef.type?.type_data === 'str') || quantityDef.type === undefined) {
    return <AutoComplete
      freeSolo
      options={suggestions || []}
      onChange={(event, value) => handleChange(value)}
      onBlur={(event) => handleChange(event.target.value)}
      ListboxProps={{style: {maxHeight: '150px'}}}
      value={value || null}
      renderInput={params => (
        <TextFieldWithHelp
          {...params}
          variant='filled' size='small' fullWidth
          {...fieldProps}
          {...otherProps}
          label={label}
        />
      )}
    />
  }

  return <TextFieldWithHelp
    select variant='filled' size='small' withOtherAdornment fullWidth
    value={value || ''}
    onChange={event => handleChange(event.target.value)}
    {...fieldProps}
    {...otherProps}
  >
    {quantityDef.type?.type_data.map(item => (
      <MenuItem value={item} key={item}>
        {item}
      </MenuItem>
    ))}
  </TextFieldWithHelp>
})
EnumEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func,
  suggestions: PropTypes.arrayOf(PropTypes.string)
}
