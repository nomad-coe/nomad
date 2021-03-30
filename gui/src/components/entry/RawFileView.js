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
import React, { useContext, useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Typography, makeStyles } from '@material-ui/core'
import { apiContext } from '../apiv1'
import { domains } from '../domains'
import { EntryPageContent } from './EntryPage'
import { errorContext } from '../errors'

const useStyles = makeStyles(theme => ({
  error: {
    marginTop: theme.spacing(2)
  }
}))

export default function RawFileView({uploadId, entryId}) {
  const classes = useStyles()
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [state, setState] = useState({entryData: null, doesNotExist: false})

  useEffect(() => {
    setState({entryData: null, doesNotExist: false})
  }, [setState, uploadId, entryId])

  useEffect(() => {
    api.entry(entryId).then(entry => {
      setState({entryData: entry.data, doesNotExist: false})
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setState({entryData: null, doesNotExist: true})
      } else {
        setState({entryData: null, doesNotExist: false})
        raiseError(error)
      }
    })
  }, [api, raiseError, entryId, setState])

  const entryData = state.entryData || {uploadId: uploadId, entryId: entryId}
  const domain = entryData.domain && domains[entryData.domain]

  if (state.doesNotExist) {
    return <EntryPageContent>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </EntryPageContent>
  }

  return (
    <EntryPageContent maxWidth={'1024px'} width={'100%'} minWidth={'800px'}>
      {domain && <domain.EntryRawView data={entryData} entryId={entryId} uploadId={uploadId} />}
    </EntryPageContent>
  )
}

RawFileView.propTypes = {
  uploadId: PropTypes.string.isRequired,
  entryId: PropTypes.string.isRequired
}
