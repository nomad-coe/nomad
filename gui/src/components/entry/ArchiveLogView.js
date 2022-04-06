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
import React, { useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { Typography, Accordion, AccordionSummary, AccordionDetails, makeStyles, FormGroup } from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import { amber } from '@material-ui/core/colors'
import { maxLogsToShow } from '../../config'
import Page from '../Page'
import { useErrors } from '../errors'
import { useApi } from '../api'
import { useEntryContext } from './EntryContext'
import Checkbox from '@material-ui/core/Checkbox'
import FormControlLabel from '@material-ui/core/FormControlLabel'

const useLogEntryStyles = makeStyles(theme => ({
  warning: {
    color: amber[700]
  },
  exception: {
    overflowX: 'scroll',
    margin: 0
  }
}))


const LogEntry = React.memo(function LogEntry(props) {
  const classes = useLogEntryStyles()
  const {entry} = props
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
    </Accordion>
  )
})
LogEntry.propTypes = {
  entry: PropTypes.object.isRequired
}

const useStyles = makeStyles(theme => ({
  moreLogs: {
    marginTop: theme.spacing(2)
  },
  downloadFab: {
    zIndex: 1,
    right: 32,
    bottom: 32,
    position: 'fixed !important'
  }
}))

export default function ArchiveLogView(props) {
  const classes = useStyles()
  const {entryId} = useEntryContext()
  const {api} = useApi()
  const {raiseError} = useErrors()

  const [data, setData] = useState(null)
  const [doesNotExist, setDoesNotExist] = useState(false)

  const [checkList, setCheckList] = useState({
    DEBUG: true,
    ERROR: true,
    CRITICAL: true,
    WARNING: true,
    INFO: true
  }); 

  useEffect(() => {
    api.post(`/entries/${entryId}/archive/query`, {required: {processing_logs: '*'}})
      .then(response => {
        const data = response.data.archive.processing_logs || []
        setData(data)
      })
      .catch(error => {
        if (error.name === 'DoesNotExist') {
          setDoesNotExist(true)
        } else {
          raiseError(error)
        }
      })
  }, [setData, setDoesNotExist, api, raiseError, entryId])

  if (doesNotExist) {
    return (
      <Page>
        <Typography>
          No archive log does exist for this entry. Most likely the entry itself does not
          exist.
        </Typography>
      </Page>
    )
  }

  let content = 'loading ...'
  if (data) {
    content = 
    <div>
      <FormGroup row>
        <Typography style={{padding: '8px', textAlign: 'bottom'}}>
          Filter by: 
        </Typography>
        <FormControlLabel 
          control={<Checkbox 
            checked={checkList.INFO} 
            onChange={(e) => setCheckList({...checkList, [e.target.name]:e.target.checked})} 
            name={'INFO'} 
            id={'1'} />
          } 
          label={'INFO'}
          />
        <FormControlLabel 
          control={<Checkbox 
            checked={checkList.WARNING} 
            onChange={(e) => setCheckList({...checkList, [e.target.name]:e.target.checked})} 
            name={'WARNING'} 
            id={'2'}/>
          } 
          label={'WARNING'}
        />
        <FormControlLabel 
          control={<Checkbox 
            checked={checkList.ERROR} 
            onChange={(e) => setCheckList({...checkList, [e.target.name]:e.target.checked})} 
            name={'ERROR'}
            id={'3'}/>
          } 
          label={'ERROR'}
        />
        <FormControlLabel 
          control={<Checkbox 
            checked={checkList.CRITICAL} 
            onChange={(e) => setCheckList({...checkList, [e.target.name]:e.target.checked})} 
            name={'CRITICAL'}
            id={'4'}/>
          } 
          label={'CRITICAL'}
        />
        <FormControlLabel 
          control={<Checkbox 
            checked={checkList.DEBUG} 
            onChange={(e) => setCheckList({...checkList, [e.target.name]:e.target.checked})} 
            name={'DEBUG'} 
            id={'5'}/>
          } 
          label={'DEBUG'}
        />
      </FormGroup>
      {data.slice(0, maxLogsToShow).map((entry, i) => (checkList[entry.level] ? <LogEntry key={i} entry={entry}/> : null))}
      {data.length > maxLogsToShow && <Typography classes={{root: classes.moreLogs}}>
        There are {data.length - maxLogsToShow} more log entries. Download the log to see all of them.
      </Typography>}
    </div>
  }

  return (
    <Page limitedWidth>
      {content}
    </Page>
  )
}
ArchiveLogView.propTypes = {}
