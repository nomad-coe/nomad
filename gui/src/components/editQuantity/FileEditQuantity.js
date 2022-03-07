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
import React, {useCallback } from 'react'
import PropTypes from 'prop-types'
import { makeStyles, IconButton, Tooltip, TextField } from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import Dropzone from 'react-dropzone'
import { useApi } from '../api'
import { ItemButton } from '../archive/Browser'
import { useEntryContext } from '../entry/EntryContext'
import { useErrors } from '../errors'

const useFileEditQuantityStyles = makeStyles(theme => ({
  dropzone: {
    border: 0
  },
  dropzoneActive: {
    backgroundColor: theme.palette.primary.main,
    borderTopLeftRadius: 4,
    borderTopRightRadius: 4,
    '& input, & button, & label, & a': {
      color: theme.palette.primary.contrastText
    }
  }
}))
const FileEditQuantity = React.memo(({onChange, quantityDef, section, ...otherProps}) => {
  const classes = useFileEditQuantityStyles()
  const {uploadId, metadata} = useEntryContext()
  const {api} = useApi()
  const {raiseError} = useErrors()

  const value = section[quantityDef.name]

  const handleDrop = useCallback(files => {
    if (!files[0]?.name) {
      return // Not dropping a file, but something else. Ignore.
    }
    const mainfilePathSegments = metadata.mainfile.split('/')
    const mainfileDir = mainfilePathSegments.slice(0, mainfilePathSegments.length - 1).join('/')
    const mainfileDirEncoded = mainfileDir.split('/').map(segment => encodeURIComponent(segment)).join('/')
    const fullPath = mainfileDir ? mainfileDir + '/' + files[0].name : files[0].name

    const formData = new FormData() // eslint-disable-line no-undef
    formData.append('file', files[0])

    api.put(
      `/uploads/${uploadId}/raw/${mainfileDirEncoded}?wait_for_processing=true`,
      formData, {
        onUploadProgress: (progressEvent) => {
          // TODO: would be nice to show progress somehow
        }
      }
    ).catch(raiseError)
    if (onChange) {
      onChange(fullPath, section, quantityDef)
    }
  }, [api, raiseError, uploadId, metadata, onChange, section, quantityDef])

  const handleChange = useCallback(event => {
    const value = event.target.value
    if (onChange) {
      onChange(value, section, quantityDef)
    }
  }, [onChange, section, quantityDef])

  console.log('####', value)
  return (
    <Dropzone
      className={classes.dropzone} activeClassName={classes.dropzoneActive}
      onDrop={handleDrop} disableClick
    >
      <TextField
        value={value || ''} onChange={handleChange}
        size="small" variant="filled" fullWidth
        label={quantityDef.name}
        {...otherProps}
        InputProps={{
          endAdornment: (
            <React.Fragment>
              <Dropzone className={classes.dropzone} onDrop={handleDrop}>
                <IconButton size="small">
                  <Tooltip title="upload file (click or drop file)">
                    <UploadIcon/>
                  </Tooltip>
                </IconButton>
              </Dropzone>
              <ItemButton size="small" itemKey={quantityDef.name} />
            </React.Fragment>
          )
        }}
      />
    </Dropzone>
  )
})
FileEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
export default FileEditQuantity
