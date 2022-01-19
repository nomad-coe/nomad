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
import { Typography, makeStyles, Card, CardHeader, CardContent } from '@material-ui/core'
import { errorContext } from '../errors'
import { useApi } from '../api'
import RawFiles from './RawFiles'
import Page from '../Page'

const useStyles = makeStyles(theme => ({
  error: {
    marginTop: theme.spacing(2)
  }
}))

export default function RawFileView({entryId}) {
  const classes = useStyles()
  const {raiseError} = useContext(errorContext)
  const [state, setState] = useState({entryData: null, doesNotExist: false})
  const {api} = useApi()

  useEffect(() => {
    setState({entryData: null, doesNotExist: false})
  }, [setState, entryId])

  useEffect(() => {
    api.get(`/entries/${entryId}`).then(entry => {
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

  const entryData = state.entryData || {entryId: entryId}

  if (state.doesNotExist) {
    return <Page>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </Page>
  }

  return (
    <Page limitedWidth>
      <Card className={classes.root}>
        <CardHeader title="Raw files" />
        <CardContent>
          <RawFiles data={entryData} entryId={entryId} />
        </CardContent>
      </Card>
    </Page>
  )
}

RawFileView.propTypes = {
  entryId: PropTypes.string.isRequired
}
