import React from 'react'
import PropTypes from 'prop-types'
import Button from '@material-ui/core/Button'
import Dialog from '@material-ui/core/Dialog'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogTitle from '@material-ui/core/DialogTitle'
import { FormGroup, Checkbox, FormLabel } from '@material-ui/core'
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
    const { onConfirm, onClose, open, title, content,confirmLabel } = this.props
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
