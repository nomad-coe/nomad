import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Fab } from '@material-ui/core'
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
    downloadFab: {
      zIndex: 1,
      right: 32,
      bottom: 32,
      position: 'fixed !important'
    }
  });

  constructor(props) {
    super(props)
    this.state = {
      data: null
    }
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api) {
      this.update()
    }
  }

  update() {
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
        <pre>{data || 'empty log'}</pre>

        <Download
          classes={{root: classes.downloadFab}} tooltip="download logfile"
          component={Fab} className={classes.downloadFab} color="primary" size="medium"
          url={`archive/logs/${uploadId}/${calcId}`} fileName={`${calcId}.log`}
        >
          <DownloadIcon />
        </Download>
      </div>
    )
  }
}

export default compose(withApi(false, true), withStyles(ArchiveLogView.styles))(ArchiveLogView)
