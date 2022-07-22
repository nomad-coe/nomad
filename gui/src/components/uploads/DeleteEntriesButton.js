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
import React, { useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import { Tooltip, IconButton, Dialog, DialogContent, DialogContentText, DialogActions, Button } from '@material-ui/core'
import DeleteIcon from '@material-ui/icons/Delete'
import { useUploadPageContext } from './UploadPageContext'
import {useApi} from '../api'
import {useErrors} from '../errors'
import {pluralize} from '../../utils'

const DeleteEntriesButton = React.memo(({tooltip, disabled, buttonProps, dark, selectedEntries, selectedCount, setSelected}) => {
  const {uploadId, updateUpload} = useUploadPageContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [openDeleteConfirmDialog, setOpenDeleteConfirmDialog] = useState(false)

  const handleClick = useCallback(() => {
    setOpenDeleteConfirmDialog(true)
  }, [setOpenDeleteConfirmDialog])

  const handleDelete = useCallback((includeParentFolders) => {
    setOpenDeleteConfirmDialog(false)
    const requestBody = {query: selectedEntries, include_parent_folders: includeParentFolders}
    api.post(`uploads/${uploadId}/action/delete-entry-files`, requestBody)
      .then(results => {
        updateUpload({upload: results.data})
        setSelected([])
      })
      .catch(err =>
        raiseError(err)
      )
  }, [setOpenDeleteConfirmDialog, api, raiseError, uploadId, selectedEntries, setSelected, updateUpload])

  return (
    <React.Fragment>
      <IconButton
        {...buttonProps}
        disabled={disabled}
        onClick={handleClick}
        style={dark ? {color: 'white'} : null}
      >
        <Tooltip title={tooltip || 'Delete selected entries'}>
          <DeleteIcon />
        </Tooltip>
      </IconButton>
      <Dialog
        open={openDeleteConfirmDialog}
        aria-describedby="alert-dialog-description"
      >
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            <b>{`Really delete selected ${pluralize('entry', selectedCount, true)}?`}</b>
          </DialogContentText>
          <DialogContentText>
            You can choose to delete only the mainfiles, or to delete the mainfiles and their folders.
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOpenDeleteConfirmDialog(false)} autoFocus>Cancel</Button>
          <Button onClick={() => handleDelete(false)}>Delete mainfiles</Button>
          <Button onClick={() => handleDelete(true)}>Delete mainfiles and folders</Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>)
})
DeleteEntriesButton.propTypes = {
  selectedEntries: PropTypes.object.isRequired,
  selectedCount: PropTypes.number.isRequired,
  setSelected: PropTypes.func.isRequired,
  tooltip: PropTypes.string,
  disabled: PropTypes.bool,
  buttonProps: PropTypes.object,
  dark: PropTypes.bool
}

export default DeleteEntriesButton
