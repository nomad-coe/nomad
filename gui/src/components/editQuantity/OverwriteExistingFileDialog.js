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
import React from 'react'
import PropTypes from 'prop-types'
import {Button, Dialog, DialogActions, DialogContent, DialogTitle, makeStyles} from '@material-ui/core'
import DialogContentText from '@material-ui/core/DialogContentText'

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%',
    minWidth: 400
  }
}))

const OverwriteExistingFileDialog = React.memo(({open, data, message, onOverwrite, onCancel}) => {
  const classes = useStyles()

  return <React.Fragment>
    <Dialog classes={{paper: classes.dialog}} open={open} disableEscapeKeyDown data-testid='overwrite-reference-dialog'>
      <DialogTitle>Overwrite file</DialogTitle>
      <DialogContent>
        <DialogContentText>
          {message || `There is already a file with the same name in this path. Do you want to overwrite the file?`}
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={onCancel} color="secondary">
          No
        </Button>
        <Button onClick={() => onOverwrite(data)} color="secondary">
          Yes
        </Button>
      </DialogActions>
    </Dialog>
  </React.Fragment>
})
OverwriteExistingFileDialog.propTypes = {
  open: PropTypes.bool,
  data: PropTypes.any,
  message: PropTypes.string,
  onOverwrite: PropTypes.func,
  onCancel: PropTypes.func
}

export default OverwriteExistingFileDialog
