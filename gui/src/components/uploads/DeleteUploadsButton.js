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
import React, {useState, useCallback, useEffect} from 'react'
import PropTypes from 'prop-types'
import { Tooltip, IconButton, Dialog, DialogContent, DialogContentText, DialogActions, Button } from '@material-ui/core'
import DeleteIcon from '@material-ui/icons/Delete'
import {pluralize} from '../../utils'
import {useApi} from '../api'
import {useErrors} from '../errors'
import DeletingReferencesTable from './DeletingReferencesTable'

const DeleteUploadsButton = React.memo(({tooltip, disabled, buttonProps, dark, selectedUploads, selectedCount, onSubmitted}) => {
  const [openDeleteConfirmDialog, setOpenDeleteConfirmDialog] = useState(false)
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [entryReferences, setEntryReferences] = useState([])
  const [brokenEntries, setBrokenEntries] = useState([])

  useEffect(() => {
    if (!openDeleteConfirmDialog) {
      return
    }
    const requestBody = {owner: 'visible', query: {entry_references: {'target_upload_id:any': selectedUploads || []}}}
    api.post(`entries/query`, requestBody)
      .then(results => {
        const data = results?.data || []
        if (data.length > 0) {
          const brokenEntries = data.filter(entry => !selectedUploads.includes(entry.upload_id))
          setBrokenEntries(brokenEntries)
          const allReferences = data.map(entry => entry.entry_references).flat()
          const references = [...new Map(allReferences.map(entry => [entry.target_entry_id, entry])).values()]
          setEntryReferences(references)
        }
      })
      .catch(err =>
        raiseError(err)
      )
  }, [api, openDeleteConfirmDialog, raiseError, selectedUploads])

  const handleClick = useCallback(() => {
    setOpenDeleteConfirmDialog(true)
  }, [setOpenDeleteConfirmDialog])

  const handleDelete = useCallback(() => {
    setOpenDeleteConfirmDialog(false)
    onSubmitted()
  }, [onSubmitted])

  return (
    <React.Fragment>
      <IconButton
        {...buttonProps}
        disabled={disabled}
        onClick={handleClick}
        style={dark ? {color: 'white'} : null}
      >
        <Tooltip title={tooltip || 'Delete selected uploads'}>
          <DeleteIcon />
        </Tooltip>
      </IconButton>
      <Dialog
        open={openDeleteConfirmDialog}
      >
        <DialogContent>
          <DialogContentText>
            <b>{`Please confirm deleting the ${pluralize('upload', selectedCount, false)}.`}</b>
          </DialogContentText>
          <DialogContentText>
            {`You have selected ${pluralize('upload', selectedCount, true)}. Are you sure you want to delete the selected ${pluralize('upload', selectedCount, false)}?`}
          </DialogContentText>
          <DeletingReferencesTable entryReferences={entryReferences} brokenEntries={brokenEntries}/>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOpenDeleteConfirmDialog(false)} autoFocus>Cancel</Button>
          <Button onClick={() => handleDelete()}>{brokenEntries.length > 0 ? 'Delete anyway' : 'Delete'}</Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>)
})
DeleteUploadsButton.propTypes = {
  selectedUploads: PropTypes.arrayOf(PropTypes.string).isRequired,
  selectedCount: PropTypes.number.isRequired,
  onSubmitted: PropTypes.func.isRequired,
  tooltip: PropTypes.string,
  disabled: PropTypes.bool,
  buttonProps: PropTypes.object,
  dark: PropTypes.bool
}

export default DeleteUploadsButton
