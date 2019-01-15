import React from 'react'
import PropTypes from 'prop-types'
import { MuiThemeProvider, withStyles } from '@material-ui/core'
import RepoIcon from '@material-ui/icons/Cloud'
import ArchiveIcon from '@material-ui/icons/Storage'
import ArchiveLogIcon from '@material-ui/icons/BugReport'
import { repoTheme, archiveTheme } from '../config'
import CalcDialog from './CalcDialog'
import RepoCalcView from './RepoCalcView'
import ArchiveCalcView from './ArchiveCalcView'
import ArchiveLogView from './ArchiveLogView'

class CalcDialogButtons extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.string,
    calcId: PropTypes.string,
    disabled: PropTypes.bool,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      overflow: 'hidden',
      whiteSpace: 'nowrap'
    }
  });

  onClickOpen

  render() {
    const { classes, ...others } = this.props

    return (
      <div className={classes.root}>
        <MuiThemeProvider theme={repoTheme}>
          <CalcDialog icon={<RepoIcon/>} component={RepoCalcView} title={'Repository calculation data'} {...others} />
        </MuiThemeProvider>
        <MuiThemeProvider theme={archiveTheme}>
          <CalcDialog icon={<ArchiveIcon/>} component={ArchiveCalcView} title={'Archive calculation data'} {...others} />
          <CalcDialog icon={<ArchiveLogIcon/>} component={ArchiveLogView} title={'Calculation processing logs'} {...others} />
        </MuiThemeProvider>
      </div>
    )
  }
}

export default withStyles(CalcDialogButtons.styles)(CalcDialogButtons)
