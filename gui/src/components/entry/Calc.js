import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Tab, Tabs } from '@material-ui/core'
import ArchiveCalcView from './ArchiveCalcView'
import ArchiveLogView from './ArchiveLogView'
import RepoCalcView from './RepoCalcView'

class Calc extends React.Component {
  static styles = theme => ({
    root: {
    },
    content: {
      padding: 0,
      maxWidth: 1024,
      margin: 'auto'
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
  }

  state = {
    viewIndex: 0
  }

  render() {
    const { classes, ...calcProps } = this.props
    const { viewIndex } = this.state

    return (
      <div className={classes.root}>
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

        <div className={classes.content}>
          <div style={viewIndex !== 0 ? {display: 'none'} : {}} >
            <RepoCalcView {...calcProps} />
          </div>
          <div style={viewIndex !== 1 ? {display: 'none'} : {}} >
            <ArchiveCalcView {...calcProps} />
          </div>
          <div style={viewIndex !== 2 ? {display: 'none'} : {}} >
            <ArchiveLogView {...calcProps} />
          </div>
        </div>
      </div>
    )
  }
}

export default withStyles(Calc.styles)(Calc)
