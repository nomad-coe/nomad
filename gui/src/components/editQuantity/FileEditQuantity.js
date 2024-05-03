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
import React, {useCallback, useRef, useState} from 'react'
import PropTypes from 'prop-types'
import {
  makeStyles, IconButton, Tooltip, TextField
} from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import { useDropzone } from 'react-dropzone'
import { useApi } from '../api'
import {ItemButton, useLane} from '../archive/Browser'
import { useEntryStore } from '../entry/EntryContext'
import OverwriteExistingFileDialog from './OverwriteExistingFileDialog'
import UploadProgressDialog from '../uploads/UploadProgressDialog'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

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
const FileEditQuantity = React.memo(props => {
  const classes = useFileEditQuantityStyles()
  const {onChange, onFailed, quantityDef, value, ...otherProps} = props
  const {index} = otherProps
  const {uploadId, metadata} = useEntryStore() || {}
  const {api} = useApi()
  const [askForOverwrite, setAskForOverwrite] = useState(false)
  const dropedFiles = useRef([])
  const [uploading, setUploading] = useState(null)
  const lane = useLane()
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const uploadFile = useCallback((files, overwrite = false) => {
    const mainfilePathSegments = metadata.mainfile.split('/')
    const mainfileDir = mainfilePathSegments.slice(0, mainfilePathSegments.length - 1).join('/')
    const mainfileDirEncoded = mainfileDir.split('/').map(segment => encodeURIComponent(segment)).join('/')
    let fullPaths = files.map((file) => mainfileDir ? mainfileDir + '/' + file.name : file.name)

    const formData = new FormData() // eslint-disable-line no-undef
    if (quantityDef.shape.length > 0) {
      files.forEach((file) => {
        formData.append('file', file)
      })
    } else {
      formData.append('file', files[0])
      fullPaths = fullPaths[0]
    }

    return new Promise((resolve, reject) => {
      api.put(
        `/uploads/${uploadId}/raw/${mainfileDirEncoded}?overwrite_if_exists=${overwrite}`,
        formData, {
          onUploadProgress: (progressEvent) => {
            const percentCompleted = Math.round((progressEvent.loaded * 100) / progressEvent.total)
            setUploading(percentCompleted)
          }
        }
      )
        .then(response => resolve({response, fullPaths}))
        .catch(error => reject(error))
        .finally(() => {
          setUploading(null)
        })
    })
  }, [api, uploadId, metadata, quantityDef])

  const handleDropFiles = useCallback(files => {
    if (!files?.[0]?.name) {
      return // Not dropping a file, but something else. Ignore.
    }

    dropedFiles.current = files
    uploadFile(files).then(({response, fullPaths}) => {
      if (onChange) {
        onChange(fullPaths)
      }
    }).catch(error => {
      if (error.apiMessage === "The provided path already exists and overwrite_if_exists is set to False.") {
        setAskForOverwrite(true)
      } else {
        onFailed(new Error(error))
      }
    })
  }, [uploadFile, onChange, onFailed])

  const handleChange = useCallback(event => {
    const value = event.target.value
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  const {getRootProps, getInputProps, open, isDragAccept} = useDropzone({onDrop: handleDropFiles, noClick: true})
  const dropzoneClassName = isDragAccept ? classes.dropzoneActive : classes.dropzone

  const handleOverwrite = useCallback((files) => {
    if (!files?.[0]?.name) {
      return
    }

    uploadFile(files, true)
      .then(({response, fullPaths}) => {
        if (onChange) {
          onChange(fullPaths)
        }
      })
      .catch(error => {
        onFailed(new Error(error))
      })
    setAskForOverwrite(false)
  }, [onChange, onFailed, uploadFile])

  return (
    <div {...getRootProps({className: dropzoneClassName})}>
      <UploadProgressDialog uploading={uploading} />
      <input {...getInputProps()} />
      <TextField
        value={value || ''} onChange={handleChange}
        size="small" variant="filled" fullWidth
        label={label}
        {...otherProps}
        InputProps={{
          endAdornment: (
            <React.Fragment>
              <IconButton size="small" onClick={open} >
                <Tooltip title="upload file (click or drop file)">
                  <UploadIcon/>
                </Tooltip>
              </IconButton>
              {lane && value && (
                <ItemButton size="small" itemKey={index === undefined ? quantityDef.name : `${quantityDef.name}/${index}`} />
              )}
            </React.Fragment>
          )
        }}
      />
      <OverwriteExistingFileDialog
        open={askForOverwrite}
        data={dropedFiles.current}
        onOverwrite={handleOverwrite}
        onCancel={() => setAskForOverwrite(false)}
      />
    </div>
  )
})
FileEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func,
  onFailed: PropTypes.func
}
export default FileEditQuantity
