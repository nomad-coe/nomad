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
import Button from '@material-ui/core/Button'
import Dialog from '@material-ui/core/Dialog'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogContentText from '@material-ui/core/DialogContentText'
import DialogTitle from '@material-ui/core/DialogTitle'

class CommingSoon extends React.Component {
  static propTypes = {
    children: PropTypes.any,
    onClose: PropTypes.func.isRequired,
    open: PropTypes.bool,
    title: PropTypes.string
  }

  render() {
    return (
      <Dialog
        open={this.props.open}
        onClose={this.handleClose}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">{this.props.title || 'Comming soon'}</DialogTitle>
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            {this.props.children || 'This feature is not yet implemented, but will be available soon.'}
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={this.props.onClose} color="primary" autoFocus>
            Got it!
          </Button>
        </DialogActions>
      </Dialog>
    )
  }
}

export default CommingSoon
