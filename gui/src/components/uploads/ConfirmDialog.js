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
            <Markdown>{`
              If you agree the selected uploads will move out of your private staging
              area into the public [NOMAD Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/).
              If you wish to put an embargo on your data it will last upto 36 month. Afterwards, your data will
              be made public. All public data will be made available under the Creative
              Commons Attribution license ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).

              The published data will be added to the NOMAD Repository's index overnight.
              Therefore, it will take until tomorrow before your data appears in the
              [NOMAD Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/).
            `}</Markdown>

            <FormGroup row style={{alignItems: 'center'}}>
              <Checkbox
                checked={!withEmbargo}
                onChange={() => this.setState({withEmbargo: !withEmbargo})}
              />
              <FormLabel>publish without embargo</FormLabel>
            </FormGroup>
          </DialogContent>
          <DialogActions>
            <Button onClick={onClose}>
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
