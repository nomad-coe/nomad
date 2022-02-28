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
import React, {useCallback, useContext } from 'react'
import PropTypes from 'prop-types'
import { makeStyles, Grid, IconButton, Tooltip } from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import Dropzone from 'react-dropzone'
import { useApi } from '../api'
import { browserContext } from '../archive/Browser'
import { useEntryContext } from '../entry/EntryContext'
import { StringEditQuantity } from './EditQuantity'

const useFileEditQuantityStyles = makeStyles(theme => ({
  dropzone: {
    border: 0
  },
  dropzoneActive: {
    backgroundColor: theme.palette.primary.light
  }
}))
const FileEditQuantity = React.memo((props) => {
  const classes = useFileEditQuantityStyles()
  const browser = useContext(browserContext)
  const {uploadId, metadata} = useEntryContext()
  const { api } = useApi()

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

    browser.blockUntilProcessed({
      uploadId: uploadId,
      apiCall: api.put(`/uploads/${uploadId}/raw/${mainfileDirEncoded}`, formData, {
        onUploadProgress: (progressEvent) => {
          // TODO: would be nice to show progress somehow
        }
      }),
      apiCallText: 'Uploading file',
      onSuccess: () => {
        if (props.onChange) {
          props.onChange(fullPath, props.section, props.quantityDef)
        }
        browser.lanes.current.forEach(lane => lane.adaptor?.onFilesUpdated(uploadId, mainfileDir))
        browser.update()
      }
    })
  }, [api, browser, uploadId, metadata, props])

  return (
    <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
      <Grid item style={{flexGrow: 1}}>
        <Dropzone
          className={classes.dropzone} activeClassName={classes.dropzoneActive}
          onDrop={handleDrop} disableClick
        >
          <StringEditQuantity {...props} />
        </Dropzone>
      </Grid>
      <Grid item style={{marginTop: 'auto', marginBottom: 'auto'}}>
        <Dropzone className={classes.dropzone} onDrop={handleDrop}>
          <IconButton size="small">
            <Tooltip title="upload file (click or drop file)">
              <UploadIcon/>
            </Tooltip>
          </IconButton>
        </Dropzone>
      </Grid>
    </Grid>)
})
FileEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
export default FileEditQuantity
