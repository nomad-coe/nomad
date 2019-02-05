import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Dialog, DialogContent, DialogActions, Button, DialogTitle, Tab, Tabs } from '@material-ui/core'
import SwipeableViews from 'react-swipeable-views'
import ArchiveCalcView from './ArchiveCalcView'
import ArchiveLogView from './ArchiveLogView'
import RepoCalcView from './RepoCalcView'

class CalcDialog extends React.Component {
  static styles = theme => ({
    dialog: {

    },
    dialogTitle: {
      padding: 0
    },
    dialogContent: {
      padding: 0
    },
    tabContent: {
      padding: `0 ${theme.spacing.unit * 3}px`,
      overflowY: 'auto',
      height: '70vh',
      zIndex: 1,
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
    open: false,
    viewIndex: 0
  }

  render() {
    const { classes, onClose, ...calcProps } = this.props
    const { viewIndex } = this.state

    return (
      <Dialog className={classes.dialog} open={true} onClose={onClose} fullWidth={true} maxWidth={'md'} >
        <DialogTitle disableTypography classes={{root: classes.dialogTitle}}>
          <Tabs
            className={classes.tabs}
            value={viewIndex}
            onChange={(event, state) => this.setState({viewIndex: state})}
            indicatorColor="primary"
            textColor="primary"
            variant="fullWidth"
          >
            <Tab label="Raw data" />
            <Tab label="Archive" />
            <Tab label="Logs" />
          </Tabs>
        </DialogTitle>
        <DialogContent classes={{root: classes.dialogContent}}>
          <SwipeableViews
            index={viewIndex}
            onChangeIndex={() => null}
          >
            <div className={classes.tabContent}>
              <RepoCalcView {...calcProps} />
            </div>
            <div className={classes.tabContent}>
              <ArchiveCalcView {...calcProps} />
            </div>
            <div className={classes.tabContent}>
              <ArchiveLogView {...calcProps} />
            </div>
          </SwipeableViews>
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

export default withStyles(CalcDialog.styles)(CalcDialog)
