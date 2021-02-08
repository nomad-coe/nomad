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
import { apiContext } from '../api'
import { domains } from '../domains'
import { EntryPageContent } from './EntryPage'
import { errorContext } from '../errors'

const useStyles = makeStyles(theme => ({
  error: {
    marginTop: theme.spacing(2)
  }
}))

export default function RawFileView({uploadId, calcId}) {
  const classes = useStyles()
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [state, setState] = useState({calcData: null, doesNotExist: false})

  useEffect(() => {
    setState({calcData: null, doesNotExist: false})
  }, [setState, uploadId, calcId])

  useEffect(() => {
    api.repo(uploadId, calcId).then(data => {
      setState({calcData: data, doesNotExist: false})
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setState({calcData: null, doesNotExist: true})
      } else {
        setState({calcData: null, doesNotExist: false})
        raiseError(error)
      }
    })
  }, [api, raiseError, uploadId, calcId, setState])

  const calcData = state.calcData || {uploadId: uploadId, calcId: calcId}
  const domain = calcData.domain && domains[calcData.domain]

  if (state.doesNotExist) {
    return <EntryPageContent>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </EntryPageContent>
  }

  return (
    <EntryPageContent maxWidth={'1024px'} width={'100%'} minWidth={'800px'}>
      {domain && <domain.EntryRawView data={calcData} calcId={calcId} uploadId={uploadId} />}
    </EntryPageContent>
  )
}

RawFileView.propTypes = {
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
