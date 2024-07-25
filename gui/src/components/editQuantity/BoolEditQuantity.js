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
  Checkbox, FormLabel, RadioGroup, Radio, FormControl
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {getFieldProps, WithHelp} from './StringEditQuantity'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const BoolEditQuantity = React.memo(function BoolEditQuantity(props) {
  const {quantityDef, value, onChange, booleanLabels, ...otherProps} = props
  const fieldProps = getFieldProps(quantityDef)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleCheckboxChange = useCallback(e => {
    const value = e.target.checked
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  const handleRadioChange = useCallback(e => {
    const value = e.target.value === 'true' ? true : (e.target.value === 'false' ? false : undefined)
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  return (
    <WithHelp {...fieldProps} label={label}>
      {booleanLabels ? (
        <FormControl>
          <FormLabel>{label}</FormLabel>
          <RadioGroup row value={String(value)} onChange={handleRadioChange}>
            {Object.keys(booleanLabels).map(key => (
              <FormControlLabel
                value={key}
                key={key}
                control={<Radio {...otherProps} />}
                label={booleanLabels[key]}
              />
            ))}
          </RadioGroup>
        </FormControl>
      ) : (
        <FormControlLabel
          control={
            <Checkbox
              onChange={handleCheckboxChange}
              color="primary"
              checked={value === true}
              indeterminate={value === undefined}
              {...otherProps}
            />
          }
          label={label}
        />
      )}
    </WithHelp>
  )
})

BoolEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.oneOf([true, false, undefined]),
  onChange: PropTypes.func,
  booleanLabels: (props, propName, componentName) => {
    const labels = props[propName]
    if (!labels) return null
    if (typeof labels !== 'object') {
      return new Error(`Invalid prop \`${propName}\` supplied to \`${componentName}\`. Expected an object.`)
    }
    const validKeys = ['true', 'false', 'undefined']
    if (!Object.keys(labels).every(key => validKeys.includes(key))) {
      return new Error(`Invalid prop \`${propName}\` supplied to \`${componentName}\`. Expected keys "true", "false", "undefined".`)
    }
    return null
  }
}
