import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Fab, Typography } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
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
    error: {
      marginTop: theme.spacing.unit * 2
    },
    downloadFab: {
      zIndex: 1,
      right: 32,
      bottom: 32,
      position: 'fixed !important'
    }
  });

  static defaultState = {
      data: null,
      doesNotExist: false
    }

  state = {...ArchiveLogView.defaultState}

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api ||
        prevProps.uploadId !== this.props.uploadId ||
        prevProps.calcId !== this.props.calcId) {
      this.setState({...ArchiveLogView.defaultState})
      this.update()
    }
  }

  update() {
    const {uploadId, calcId, api, raiseError} = this.props
    api.calcProcLog(uploadId, calcId).then(data => {
      this.setState({data: data})
    }).catch(error => {
      this.setState({data: null})
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true})
      } else {
        raiseError(error)
      }
    })
  }

  render() {
    const { classes, uploadId, calcId } = this.props
    const { data, doesNotExist } = this.state

    if (doesNotExist) {
      return (
        <Typography className={classes.error}>
          No archive log does exist for this entry. Most likely the entry itself does not
          exist.
        </Typography>
      )
    }

    return (
      <div className={classes.root}>
        <pre>{data || 'loading ...'}</pre>

        <Download
          classes={{root: classes.downloadFab}} tooltip="download logfile"
          component={Fab} className={classes.downloadFab} color="secondary" size="medium"
          url={`archive/logs/${uploadId}/${calcId}`} fileName={`${calcId}.log`}
        >
          <DownloadIcon />
        </Download>
      </div>
    )
  }
}

export default compose(withApi(false, true), withStyles(ArchiveLogView.styles))(ArchiveLogView)
