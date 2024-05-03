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
import PropTypes from 'prop-types'
import AutoComplete from '@material-ui/lab/Autocomplete'
import {getFieldProps, TextFieldWithHelp} from './StringEditQuantity'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const AutocompleteEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, ...otherProps} = props
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleChange = useCallback((event, value) => {
    value = value || event.target.value
    if (onChange) {
      onChange((value === null || value === '') ? undefined : value)
    }
  }, [onChange])

  return <AutoComplete
    options={quantityDef.type.type_data}
    onChange={handleChange}
    ListboxProps={{style: {maxHeight: '150px'}}}
    value={value || null}
    renderInput={params => (
      <TextFieldWithHelp
        {...params}
        variant='filled' size='small' fullWidth
        {...getFieldProps(quantityDef)}
        {...otherProps}
        label={label}
      />
    )}
  />
})
AutocompleteEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}
