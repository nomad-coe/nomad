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
import DialogTitle from '@material-ui/core/DialogTitle'
import Markdown from '../Markdown'

class ConfirmDialog extends React.Component {
  static propTypes = {
    title: PropTypes.string.isRequired,
    confirmLabel: PropTypes.string,
    content: PropTypes.string.isRequired,
    onConfirm: PropTypes.func.isRequired,
    onClose: PropTypes.func.isRequired,
    open: PropTypes.bool.isRequired
  }

  render() {
    const { onConfirm, onClose, open, title, content, confirmLabel } = this.props
    return (
      <div>
        <Dialog
          open={open}
          onClose={onClose}
        >
          <DialogTitle>{title}</DialogTitle>
          <DialogContent>
            <Markdown>{content}</Markdown>
          </DialogContent>
          <DialogActions>
            <Button onClick={onClose}>
              Cancel
            </Button>
            <Button onClick={onConfirm} color="primary" autoFocus>
              {confirmLabel || 'Confirm'}
            </Button>
          </DialogActions>
        </Dialog>
      </div>
    )
  }
}

export default ConfirmDialog
