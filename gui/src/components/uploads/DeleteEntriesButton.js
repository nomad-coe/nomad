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
import React, { useState, useCallback, useEffect } from 'react'
import PropTypes from 'prop-types'
import {
  Tooltip, IconButton, Dialog, DialogContent, DialogContentText, DialogActions, Button,
  Checkbox, FormControlLabel
} from '@material-ui/core'
import DeleteIcon from '@material-ui/icons/Delete'
import { useUploadPageContext } from './UploadPageContext'
import {useApi} from '../api'
import {useErrors} from '../errors'
import {pluralize} from '../../utils'
import {useEntryStore} from '../entry/EntryContext'
import {useHistory} from 'react-router-dom'
import DeletingReferencesTable from './DeletingReferencesTable'

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
  const [entryReferences, setEntryReferences] = useState([])
  const [brokenEntries, setBrokenEntries] = useState([])

  const handleClick = useCallback(() => {
    setOpenDeleteConfirmDialog(true)
  }, [setOpenDeleteConfirmDialog])

  useEffect(() => {
    if (!openDeleteConfirmDialog) {
      return
    }
    const requestBody = selectedEntries?.entry_id
        ? {owner: 'visible', query: {entry_references: {'target_entry_id:any': selectedEntries?.entry_id || []}}}
        : {owner: 'visible', query: {entry_references: {'target_upload_id:any': selectedEntries?.upload_id ? [selectedEntries?.upload_id] : []}}}
    api.post(`entries/query`, requestBody)
      .then(results => {
        const data = results?.data || []
        if (data.length > 0) {
          const brokenEntries = selectedEntries?.entry_id
              ? data.filter(entry => !selectedEntries?.entry_id.includes(entry.entry_id))
              : data.filter(entry => entry.upload_id !== selectedEntries?.upload_id)
          setBrokenEntries(brokenEntries)
          const allReferences = data.map(entry => entry.entry_references).flat()
          const references = [...new Map(allReferences.map(entry => [entry.target_entry_id, entry])).values()]
          setEntryReferences(references)
        }
      })
      .catch(err =>
        raiseError(err)
      )
  }, [api, openDeleteConfirmDialog, raiseError, selectedEntries])

  const handleDelete = useCallback(() => {
    setOpenDeleteConfirmDialog(false)
    const requestBody = {query: selectedEntries || {entry_id: entryId}, include_parent_folders: deleteFolders}
    api.post(`uploads/${uploadId}/action/delete-entry-files`, requestBody)
      .then(results => {
        updateUpload && updateUpload({upload: results.data})
        setSelected && setSelected(new Set())
        entryId && history.push(`/user/uploads/upload/id/${uploadId}`)
      })
      .catch(err =>
        raiseError(err)
      )
  }, [setOpenDeleteConfirmDialog, api, raiseError, uploadId, selectedEntries, setSelected, updateUpload, entryId, history, deleteFolders])

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
          <DeletingReferencesTable entryReferences={entryReferences} brokenEntries={brokenEntries}/>
        </DialogContent>
        <DialogActions >
          <Button onClick={() => setOpenDeleteConfirmDialog(false)} autoFocus>Cancel</Button>
          <Button onClick={() => handleDelete()} data-testid={'delete-dialog-delete-button'}>
            {deleteFolders
            ? `${brokenEntries.length > 0 ? 'Delete anyway' : 'Delete'} ${pluralize('entry', count, true)} and the ${pluralize('folder', count, false)}`
            : `${brokenEntries.length > 0 ? 'Delete anyway' : 'Delete'} ${pluralize('entry', count, true)}`}
          </Button>
        </DialogActions>
      </Dialog>
    </React.Fragment> : null)
})
DeleteEntriesButton.propTypes = {
  selectedEntries: PropTypes.object,
  selectedCount: PropTypes.number,
  setSelected: PropTypes.func,
  tooltip: PropTypes.string,
  disabled: PropTypes.bool,
  buttonProps: PropTypes.object,
  dark: PropTypes.bool
}

export default DeleteEntriesButton
