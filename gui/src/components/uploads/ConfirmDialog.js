import React from 'react'
import PropTypes from 'prop-types'
import Button from '@material-ui/core/Button'
import Dialog from '@material-ui/core/Dialog'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogContentText from '@material-ui/core/DialogContentText'
import DialogTitle from '@material-ui/core/DialogTitle'
import { FormGroup, Checkbox, FormLabel } from '@material-ui/core'

class ConfirmDialog extends React.Component {
  static propTypes = {
    onPublish: PropTypes.func.isRequired,
    onClose: PropTypes.func.isRequired,
    open: PropTypes.bool.isRequired
  }

  state = {
    withEmbargo: false
  }

  render() {
    const { onPublish, onClose, open } = this.props
    const { withEmbargo } = this.state
    return (
      <div>
        <Dialog
          open={open}
          onClose={onClose}
        >
          <DialogTitle>Publish data</DialogTitle>
          <DialogContent>
            <DialogContentText>
              If you agree the selected uploads will move out of your private staging
               area into the public NOMAD.
            </DialogContentText>

            <FormGroup row style={{alignItems: 'center'}}>
              <Checkbox
                checked={!withEmbargo}
                onChange={() => this.setState({withEmbargo: !withEmbargo})}
              />
              <FormLabel>publish without embargo</FormLabel>
            </FormGroup>
          </DialogContent>
          <DialogActions>
            <Button onClick={onClose} color="primary">
              Cancel
            </Button>
            <Button onClick={() => onPublish(withEmbargo)} color="primary" autoFocus>
              {withEmbargo ? 'Publish with embargo' : 'Publish'}
            </Button>
          </DialogActions>
        </Dialog>
      </div>
    )
  }
}

export default ConfirmDialog
