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
import { domainComponents } from '../domainComponents'
import { compose } from 'recompose'
import { withApiV1 } from '../apiV1'
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
function OverviewView({uploadId, entryId, apiV1}) {
  const classes = useStyles()
  const {raiseError} = useContext(errorContext)
  const [entry, setEntry] = useState(null)
  const [exists, setExists] = useState(true)

  // When loaded for the first time, download calc data from the ElasticSearch
  // index. It is used to decide the subview to show.
  useEffect(() => {
    apiV1.entry(entryId).then(data => {
      setEntry(data)
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [apiV1, raiseError, entryId, setEntry, setExists])

  // The entry does not exist
  if (!exists) {
    return <EntryPageContent>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </EntryPageContent>
  }

  // When repo data is loaded, return a subview that depends on the domain.
  if (entry?.data?.domain) {
    const domain = entry.data.domain && domainComponents[entry.data.domain]
    return <EntryPageContent fixed>
      <domain.EntryOverview data={entry.data}/>
    </EntryPageContent>
  }
  return null
}

OverviewView.propTypes = {
  uploadId: PropTypes.string.isRequired,
  entryId: PropTypes.string.isRequired,
  apiV1: PropTypes.object.isRequired
}

export default compose(
  withApiV1(false, true)
)(OverviewView)
