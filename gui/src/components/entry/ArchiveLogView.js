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
import { Typography, Accordion, AccordionSummary, AccordionDetails, makeStyles, FormGroup, Button, Grid, FormControl, InputLabel, Input, Select, MenuItem, Chip } from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import { amber } from '@material-ui/core/colors'
import Page from '../Page'
import { useErrors } from '../errors'
import { useApi } from '../api'
import { useEntryContext } from './EntryContext'
import Checkbox from '@material-ui/core/Checkbox'
import FormControlLabel from '@material-ui/core/FormControlLabel'

const logsDefaultValues = {
  defaultLogsToShowOnFirstMount: 5,
  defaultLogsToShowOnEachMount: 10
}
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
  const {keyName} = props
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
        <Typography {...summaryProps}>{data.level}: {keyName.map((key) => `${data[key]} | `)}</Typography>
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
  entry: PropTypes.object.isRequired,
  keyName: PropTypes.array.isRequired
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
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120,
    maxWidth: 350
  },
  chips: {
    display: 'flex',
    flexWrap: 'wrap'
  },
  chip: {
    margin: 2
  },
  noLabel: {
    marginTop: theme.spacing(3)
  },
  seeMore: {
    margin: theme.spacing(3)
  }
}))

const FilterLogsByLevel = React.memo(function FilterLogsByLevel(props) {
  const {checkList} = props
  const {onCheckListChanged} = props
  return (
    <FormGroup row>
      <Typography style={{padding: '8px', textAlign: 'bottom'}}>
          Filter Logs by Level:
      </Typography>
      {Object.keys(checkList).map((key, i) => {
        return (
          <FormControlLabel key={i}
            control={<Checkbox
              checked={checkList[key]}
              onChange={onCheckListChanged}
              name={key} id={`${i}`}/>}
            label={key}
          />
        )
      })}
    </FormGroup>
  )
})
FilterLogsByLevel.propTypes = {
  checkList: PropTypes.object.isRequired,
  handleCheckListChanged: PropTypes.func.isRequired,
  onCheckListChanged: PropTypes.func.isRequired
}

const FilterLogTagsByKeys = React.memo(function FilterLogTagsByKeys(props) {
  const {className, keyName, onKeyNameChanged, uniquekeys} = props
  return (
    <FormControl className={className.formControl}>
      <InputLabel id="mutiple-chip-label">Filter keys by:</InputLabel>
      <Select
        labelId="mutiple-chip-label"
        id="mutiple-chip"
        multiple
        value={keyName}
        onChange={onKeyNameChanged}
        input={<Input id="select-multiple-chip" />}
        renderValue={(selected) => (
          <div className={className.chips}>
            {selected.map((value) => (
              <Chip key={value} label={value} className={className.chip} />
            ))}
          </div>
        )}
      >
        {uniquekeys.map((name) => (
          <MenuItem key={name} value={name}>
            {name}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})

FilterLogTagsByKeys.propTypes = {
  className: PropTypes.object.isRequired,
  keyName: PropTypes.array.isRequired,
  onKeyNameChanged: PropTypes.func.isRequired,
  uniquekeys: PropTypes.array.isRequired
}
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
  })

  const [numberOfLogs, setNumberOflogs] = useState(logsDefaultValues.defaultLogsToShowOnFirstMount)
  const [keyName, setKeyName] = useState(['parser'])
  const handleKeyNameChanged = (e) => {
    setKeyName(e.target.value)
  }

  const handleCheckListChanged = (e) => {
    setCheckList({...checkList, [e.target.name]: e.target.checked})
  }

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

    setNumberOflogs(logsDefaultValues.defaultLogsToShowOnEachMount)
  }, [setData, setDoesNotExist, api, raiseError, entryId, checkList])

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
    let uniquekeys = [...new Set(
      data.reduce((uniqueKeys, item) => [...uniqueKeys, ...Object.keys(item)], []))]
    content =
    <Grid container alignItems='center'>
      <Grid item xs={8}>
        <FilterLogsByLevel checkList={checkList} onCheckListChanged={handleCheckListChanged}/>
      </Grid>
      <Grid item xs={4} >
        <FilterLogTagsByKeys
          uniquekeys={uniquekeys}
          className={classes}
          keyName={keyName}
          onKeyNameChanged={handleKeyNameChanged}
        />
      </Grid>
      <Grid container spacing={1}>
        {data.slice(0, numberOfLogs).map((entry, i) => (checkList[entry.level]
          ? <Grid item xs={12} key={i}><LogEntry key={i} entry={entry} keyName={keyName}/></Grid> : null))}
      </Grid>
      <Grid container alignItems='center' justifyContent='center'>
        {numberOfLogs < data.length ? (<Button className={classes.seeMore} variant='contained' color='primary' onClick={() => setNumberOflogs(numberOfLogs + numberOfLogs)}>
          See More
        </Button>) : ''}
      </Grid>
    </Grid>
  }

  return (
    <Page limitedWidth>
      {content}
    </Page>
  )
}
ArchiveLogView.propTypes = {}
