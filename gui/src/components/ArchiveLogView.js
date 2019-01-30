import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, LinearProgress, Fab } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { withApi } from './api'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'

class ArchiveLogView extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {
      '& pre': {
        overflowX: 'auto'
      }
    },
    downloadFab: {
      position: 'absolute',
      zIndex: 1,
      top: theme.spacing.unit,
      right: theme.spacing.unit * 3
    }
  });

  constructor(props) {
    super(props)
    this.state = {
      data: null
    }
  }

  componentDidMount() {
    const {uploadId, calcId, api, raiseError} = this.props
    api.calcProcLog(uploadId, calcId).then(data => {
      this.setState({data: data})
    }).catch(error => {
      this.setState({data: null})
      raiseError(error)
    })
  }

  render() {
    const { classes, uploadId, calcId } = this.props
    const { data } = this.state

    return (
      <div className={classes.root}>
        <Download
          classes={{root: classes.downloadFab}} tooltip="download logfile"
          component={Fab} className={classes.downloadFab} color="primary" size="medium"
          url={`archive/logs/${uploadId}/${calcId}`} fileName={`${calcId}.log`}
        >
          <DownloadIcon />
        </Download>
        {
          data
            ? <pre>{data}</pre>
            : <LinearProgress variant="query" />
        }
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(ArchiveLogView.styles))(ArchiveLogView)
