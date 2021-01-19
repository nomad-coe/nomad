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
import { apiContext } from '../api'
import { domains } from '../domains'
import { EntryPageContent } from './EntryPage'
import { errorContext } from '../errors'
import { Typography, makeStyles } from '@material-ui/core'

const useStyles = makeStyles(theme => ({
  error: {
    marginTop: theme.spacing(2)
  },
  cardContent: {
    paddingTop: 0
  },
  topCard: {
    height: '32rem'
  },
  toggle: {
    marginBottom: theme.spacing(1)
  },
  structure: {
    marginTop: theme.spacing(1),
    width: '100%',
    height: '20rem'
  }
}))

/**
 * Shows an informative overview about the selected entry.
 */
export default function OverviewView({uploadId, calcId}) {
  const classes = useStyles()
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [repo, setRepo] = useState(null)
  const [exists, setExists] = useState(true)

  // When loaded for the first time, download calc data from the ElasticSearch
  // index. It is used quick to fetch and will be used to decide the subview to
  // show.
  useEffect(() => {
    api.repo(uploadId, calcId).then(data => {
      setRepo(data)
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, uploadId, calcId, setRepo, setExists])

  // The entry does not exist
  if (!exists) {
    return <EntryPageContent>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </EntryPageContent>
  }

  // When repo data is loaded, return a subview that depends on the domain
  if (repo) {
    const domain = repo.domain && domains[repo.domain]
    return <EntryPageContent fixed>
      {domain && <domain.EntryOverview repo={repo} uploadId={uploadId} calcId={calcId}/>}
    </EntryPageContent>
  }
  return null
}

OverviewView.propTypes = {
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
