import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Fab, Typography, ExpansionPanel, ExpansionPanelSummary, ExpansionPanelDetails } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import { amber } from '@material-ui/core/colors'
import { maxLogsToShow } from '../../config'

class LogEntryUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    entry: PropTypes.object.isRequired
  }

  static styles = theme => ({
    warning: {
      color: amber[700]
    },
    exception: {
      overflowX: 'scroll',
      margin: 0
    }
  })

  render() {
    const { classes, entry } = this.props
    const data = entry

    const summaryProps = {}
    if (data.level === 'ERROR' || data.level === 'CRITICAL') {
      summaryProps.color = 'error'
    } else if (data.level === 'WARNING') {
      summaryProps.classes = {root: classes.warning}
    }
    return (
      <ExpansionPanel>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />}>
          <Typography {...summaryProps}>{data.level}: {data.event} {(data.parser || data.normalizer) ? `(${data.parser || data.normalizer})` : ''}</Typography>
        </ExpansionPanelSummary>
        <ExpansionPanelDetails>
          <ReactJson
            src={data}
            enableClipboard={false}
            displayObjectSize={false} />
        </ExpansionPanelDetails>
        {data.exception && <ExpansionPanelDetails>
          <pre className={classes.exception}>{data.exception}</pre>
        </ExpansionPanelDetails>}
      </ExpansionPanel>)
  }
}

const LogEntry = withStyles(LogEntryUnstyled.styles)(LogEntryUnstyled)

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
      marginTop: theme.spacing(2)
    },
    moreLogs: {
      marginTop: theme.spacing(2)
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
        <Typography>
          No archive log does exist for this entry. Most likely the entry itself does not
          exist.
        </Typography>
      )
    }

    let content = 'loading ...'
    if (data) {
      content = <div>
        {data.slice(0, maxLogsToShow).map((entry, i) => <LogEntry key={i} entry={entry}/>)}
        {data.length > maxLogsToShow && <Typography classes={{root: classes.moreLogs}}>
          There are {data.length - maxLogsToShow} more log entries. Download the log to see all of them.
        </Typography>}
      </div>
    }

    return (
      <div className={classes.root}>
        {content}

        <Download
          classes={{root: classes.downloadFab}} tooltip="download logfile"
          component={Fab} className={classes.downloadFab} size="medium"
          color="primary"
          url={`archive/logs/${uploadId}/${calcId}`} fileName={`${calcId}.log`}
        >
          <DownloadIcon />
        </Download>
      </div>
    )
  }
}

export default compose(withApi(false, true), withStyles(ArchiveLogView.styles))(ArchiveLogView)
