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
  Button,
  Tooltip
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {getFieldProps} from './StringEditQuantity'
import {getDisplayLabel} from "../../utils"
import {useErrors } from '../errors'
import {useEntryStore } from '../entry/EntryContext'
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const ActionEditQuantity = React.memo(function ButtonEditQuantity(props) {
  const {quantityDef, value, onChange, ...otherProps} = props
  const fieldProps = getFieldProps(quantityDef)
  const {saveArchive} = useEntryStore()
  const {raiseError} = useErrors()
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleClick = useCallback(() => {
    onChange?.(true)
    saveArchive().catch(raiseError)
  }, [onChange, saveArchive, raiseError])

  return <Tooltip title={fieldProps.helpDescription}>
    <Button
      variant="outlined"
      disabled={value}
      onClick={handleClick}
      color="primary"
      {...otherProps}
    >{label}</Button>
  </Tooltip>
})
ActionEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.bool,
  onChange: PropTypes.func
}
