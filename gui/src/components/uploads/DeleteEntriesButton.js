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
import {
  Tooltip, IconButton, Dialog, DialogContent, DialogContentText, DialogActions, Button, Checkbox, FormControlLabel
} from '@material-ui/core'
import DeleteIcon from '@material-ui/icons/Delete'
import { useUploadPageContext } from './UploadPageContext'
import {useApi} from '../api'
import {useErrors} from '../errors'
import {pluralize} from '../../utils'
import {useEntryStore} from '../entry/EntryContext'
import {useHistory} from 'react-router-dom'

export const DeleteEntriesButton = React.memo(({tooltip, disabled, buttonProps, dark, selectedEntries, selectedCount, setSelected}) => {
  const history = useHistory()
  const uploadContext = useUploadPageContext()
  const entryContext = useEntryStore()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [openDeleteConfirmDialog, setOpenDeleteConfirmDialog] = useState(false)
  const [deleteFolders, setDeleteFolders] = useState(false)
  const {uploadId} = uploadContext || entryContext
  const {updateUpload} = uploadContext || {}
  const {entryId, editable} = entryContext || {}
  const count = entryId ? 1 : selectedCount

  const handleClick = useCallback(() => {
    setOpenDeleteConfirmDialog(true)
  }, [setOpenDeleteConfirmDialog])

  const handleDelete = useCallback((includeParentFolders) => {
    setOpenDeleteConfirmDialog(false)
    const requestBody = {query: selectedEntries || {entry_id: entryId}, include_parent_folders: includeParentFolders}
    api.post(`uploads/${uploadId}/action/delete-entry-files`, requestBody)
      .then(results => {
        updateUpload && updateUpload({upload: results.data})
        setSelected && setSelected(new Set())
        entryId && history.push(`/user/uploads/upload/id/${uploadId}`)
      })
      .catch(err =>
        raiseError(err)
      )
  }, [setOpenDeleteConfirmDialog, api, raiseError, uploadId, selectedEntries, setSelected, updateUpload, entryId, history])

  return (
    uploadContext || (entryContext && editable) ? <React.Fragment>
      <IconButton
        {...buttonProps}
        disabled={disabled}
        onClick={handleClick}
        style={dark ? {color: 'white'} : null}
      >
        <Tooltip title={tooltip || entryContext ? 'Delete entry' : 'Delete selected entries'}>
          <DeleteIcon />
        </Tooltip>
      </IconButton>
      <Dialog
        open={openDeleteConfirmDialog}
        aria-describedby="alert-dialog-description"
      >
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            <b>{`Please confirm deleting the ${pluralize('entry', count, false)}.`}</b>
          </DialogContentText>
          <DialogContentText>
            You have selected {pluralize('entry', count, true)}. Are you sure you want to delete the selected {pluralize('entry', count, false)}?
          </DialogContentText>
          <FormControlLabel
            control={<Checkbox
              onChange={event => setDeleteFolders(event.target.checked)}
              color='primary'
              checked={deleteFolders}
              name='Delete folders'
            />}
            label={`Delete also the ${pluralize("entry's", count, false)} folders. There might be more entries in them.`}
          />
        </DialogContent>
        <DialogActions >
          <Button onClick={() => setOpenDeleteConfirmDialog(false)} autoFocus>Cancel</Button>
          {!deleteFolders && <Button onClick={() => handleDelete(false)}>Delete {pluralize('entry', count, true)}</Button>}
          {deleteFolders && <Button onClick={() => handleDelete(true)}>Delete {pluralize('entry', count, true)} and the {pluralize('folder', count, false)}</Button>}
        </DialogActions>
      </Dialog>
    </React.Fragment> : null)
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
