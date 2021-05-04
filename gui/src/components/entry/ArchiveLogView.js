/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Fab, Typography, Accordion, AccordionSummary, AccordionDetails } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import Download from './Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import { amber } from '@material-ui/core/colors'
import { maxLogsToShow } from '../../config'
import { EntryPageContent } from './EntryPage'

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
      <Accordion>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography {...summaryProps}>{data.level}: {data.event} {(data.parser || data.normalizer) ? `(${data.parser || data.normalizer})` : ''}</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <ReactJson
            src={data}
            enableClipboard={false}
            displayObjectSize={false} />
        </AccordionDetails>
        {data.exception && <AccordionDetails>
          <pre className={classes.exception}>{data.exception}</pre>
        </AccordionDetails>}
      </Accordion>)
  }
}

const LogEntry = withStyles(LogEntryUnstyled.styles)(LogEntryUnstyled)

class ArchiveLogView extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    entryId: PropTypes.string.isRequired
  }

  static styles = theme => ({
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
        prevProps.entryId !== this.props.entryId) {
      this.setState({...ArchiveLogView.defaultState})
      this.update()
    }
  }

  update() {
    const {uploadId, entryId, api, raiseError} = this.props
    api.calcProcLog(uploadId, entryId).then(data => {
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
    const { classes, uploadId, entryId } = this.props
    const { data, doesNotExist } = this.state

    if (doesNotExist) {
      return (
        <EntryPageContent>
          <Typography>
            No archive log does exist for this entry. Most likely the entry itself does not
            exist.
          </Typography>
        </EntryPageContent>
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
      <EntryPageContent maxWidth={'1024px'} width={'100%'} minWidth={'800px'}>
        {content}
        <Download
          classes={{root: classes.downloadFab}} tooltip="download logfile"
          component={Fab} className={classes.downloadFab} size="medium"
          color="primary"
          url={`archive/logs/${uploadId}/${entryId}`} fileName={`${entryId}.log`}
        >
          <DownloadIcon />
        </Download>
      </EntryPageContent>
    )
  }
}

export default compose(
  withApi(false, true),
  withStyles(ArchiveLogView.styles)
)(ArchiveLogView)
