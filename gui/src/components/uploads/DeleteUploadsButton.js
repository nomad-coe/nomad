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
import React, {useState, useCallback} from 'react'
import PropTypes from 'prop-types'
import { Tooltip, IconButton, Dialog, DialogContent, DialogContentText, DialogActions, Button } from '@material-ui/core'
import DeleteIcon from '@material-ui/icons/Delete'
import {pluralize} from '../../utils'

const DeleteUploadsButton = React.memo(({tooltip, disabled, buttonProps, dark, selectedCount, onSubmitted}) => {
  const [openDeleteConfirmDialog, setOpenDeleteConfirmDialog] = useState(false)

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
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOpenDeleteConfirmDialog(false)} autoFocus>Cancel</Button>
          <Button onClick={() => handleDelete()}>Delete</Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>)
})
DeleteUploadsButton.propTypes = {
  selectedCount: PropTypes.number.isRequired,
  onSubmitted: PropTypes.func.isRequired,
  tooltip: PropTypes.string,
  disabled: PropTypes.bool,
  buttonProps: PropTypes.object,
  dark: PropTypes.bool
}

export default DeleteUploadsButton
