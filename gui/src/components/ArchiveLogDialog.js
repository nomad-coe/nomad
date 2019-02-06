import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Dialog, DialogContent, DialogActions, Button, DialogTitle } from '@material-ui/core'
import ArchiveLogView from './ArchiveLogView'

class ArchiveLogDialog extends React.Component {
  static styles = theme => ({
    dialog: {},
    dialogTitle: {},
    dialogContent: {
      position: 'relative'
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    onClose: PropTypes.func.isRequired
  }

  state = {
    open: false
  }

  render() {
    const { classes, onClose, ...calcProps } = this.props

    return (
      <Dialog className={classes.dialog} open={true} onClose={onClose} fullWidth={true} maxWidth={'md'} >
        <DialogTitle disableTypography classes={{root: classes.dialogTitle}}>
         Processing logs
        </DialogTitle>
        <DialogContent classes={{root: classes.dialogContent}}>
          <ArchiveLogView {...calcProps} />
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose} color="primary" autoFocus>
              Close
          </Button>
        </DialogActions>
      </Dialog>
    )
  }
}

export default withStyles(ArchiveLogDialog.styles)(ArchiveLogDialog)
