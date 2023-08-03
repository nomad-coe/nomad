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
import React, {useState, useCallback, useEffect, useMemo} from 'react'
import PropTypes from 'prop-types'
import {
  Tooltip,
  IconButton,
  Dialog,
  DialogContent,
  DialogContentText,
  DialogActions,
  Button
} from '@material-ui/core'
import DeleteIcon from '@material-ui/icons/Delete'
import {pluralize} from '../../utils'
import {useApi} from '../api'
import {useErrors} from '../errors'
import DeletingReferencesTable from './DeletingReferencesTable'

/**
 * Button for deleting uploads. Shows a confirmation dialog with additional
 * information about the effects of the deletion.
 */
const DeleteUploadsButton = React.memo(({
  tooltip,
  disabled,
  buttonProps,
  dark,
  uploads,
  onConfirm,
  'data-testid': testID
}) => {
  const [openDialog, setOpenDialog] = useState(false)
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [entryReferences, setEntryReferences] = useState([])
  const [brokenEntries, setBrokenEntries] = useState([])
  const {uploadIDs, nUploads, shared, uploadText, haveText} = useMemo(() => {
    const nUploads = uploads?.length || 0
    const uploadText = pluralize('upload', nUploads, false)
    const haveText = pluralize('has', nUploads, false)
    return {
      uploadIDs: uploads?.map(upload => upload.upload_id) || [],
      nUploads,
      shared: uploads?.some(entry => entry?.coauthors?.length > 0 || entry?.reviewers?.length > 0),
      uploadText,
      haveText
    }
  }, [uploads])

  // Fetch information about whether deleting the uploads will result in broken
  // references.
  useEffect(() => {
    if (!openDialog) return

    api.post(`entries/query`,
      {
        owner: 'visible',
        query: {entry_references: {'target_upload_id:any': uploadIDs}}
      }
    )
      .then(results => {
        const data = results?.data || []
        if (data.length > 0) {
          const brokenEntries = data.filter(entry => !uploadIDs.includes(entry.upload_id))
          setBrokenEntries(brokenEntries)
          const allReferences = data.map(entry => entry.entry_references).flat()
          const references = [...new Map(allReferences.map(entry => [entry.target_entry_id, entry])).values()]
          setEntryReferences(references)
        }
      })
      .catch(raiseError)
  }, [api, openDialog, raiseError, uploadIDs])

  const handleClick = useCallback(() => {
    setOpenDialog(true)
  }, [])

  const handleDelete = useCallback(() => {
    setOpenDialog(false)
    onConfirm && onConfirm()
  }, [onConfirm])

  return (
    <React.Fragment>
      <IconButton
        {...buttonProps}
        disabled={disabled}
        onClick={handleClick}
        style={dark ? {color: 'white'} : null}
        data-testid={testID}
      >
        <Tooltip title={tooltip || `Delete selected ${uploadText}`}>
          <DeleteIcon />
        </Tooltip>
      </IconButton>
      <Dialog open={openDialog}>
        <DialogContent>
          <DialogContentText>
            <b>{`Please confirm deleting the ${uploadText}`}</b>
          </DialogContentText>
          <DialogContentText>
            {nUploads > 1 ? `You have selected ${nUploads} ${uploadText}. ` : ''}
            {shared ? `The ${uploadText} you are about to delete ${haveText} been shared with at least one coauthor or reviewer. ` : ''}
            {`Are you sure you want to delete the selected ${uploadText}?`}
          </DialogContentText>
          <DeletingReferencesTable entryReferences={entryReferences} brokenEntries={brokenEntries}/>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOpenDialog(false)} autoFocus>Cancel</Button>
          <Button onClick={() => handleDelete()}>{brokenEntries.length > 0 ? 'Delete anyway' : 'Delete'}</Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>)
})
DeleteUploadsButton.propTypes = {
  uploads: PropTypes.arrayOf(PropTypes.object).isRequired,
  onConfirm: PropTypes.func.isRequired,
  tooltip: PropTypes.string,
  disabled: PropTypes.bool,
  buttonProps: PropTypes.object,
  dark: PropTypes.bool,
  'data-testid': PropTypes.string
}

export default DeleteUploadsButton
