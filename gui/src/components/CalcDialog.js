import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Dialog, DialogContent, DialogTitle, DialogActions, Button, IconButton, Typography } from '@material-ui/core'

class CalcDialog extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    disabled: PropTypes.bool,
    icon: PropTypes.element,
    component: PropTypes.func.isRequired,
    title: PropTypes.string
  }

  state = {
    open: false
  }

  render() {
    const {title, icon, component, disabled, uploadId, calcId, raiseError} = this.props
    return (
      <span>
        <IconButton color="primary" onClick={() => this.setState({open: true})} disabled={disabled}>{icon}</IconButton>
        <Dialog open={this.state.open} onClose={() => this.setState({open: false})} fullWidth={true} maxWidth={'md'}>
          <DialogTitle>
            <Typography variant="h6" gutterBottom>{title || 'Calculation data'}</Typography>
            <Typography variant="caption" gutterBottom>{`${uploadId}/${calcId}`}</Typography>
          </DialogTitle>
          <DialogContent>
            {component({uploadId: uploadId, calcId: calcId, raiseError: raiseError})}
          </DialogContent>
          <DialogActions>
            <Button onClick={() => this.setState({open: false})} color="primary" autoFocus>
                Close
            </Button>
          </DialogActions>
        </Dialog>
      </span>
    )
  }
}

export default withStyles(CalcDialog.styles)(CalcDialog)
