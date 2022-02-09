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
import React, { useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { Card, CardContent, Typography, makeStyles } from '@material-ui/core'
import ArchiveBrowser from '../archive/ArchiveBrowser'
import Page from '../Page'
import { useErrors } from '../errors'
import { useApi } from '../api'
import { ApiDataContext } from '../buttons/SourceDialogButton'

export const help = `
The NOMAD **archive** provides data and meta-data in a common hierarchical format based on
well-defined quantity definitions that we call *metainfo*. This representation
is independent from the raw data format and provides a homogenous data stock.

You can click the various quantity values to see the quantity definition. Similarly,
you can click section names to get more information. Browse the *metainfo* to
learn more about NOMAD's archive format [here](/metainfo).
`

const useStyles = makeStyles(theme => ({
  archiveBrowser: {
    marginTop: theme.spacing(2)
  },
  error: {
    marginTop: theme.spacing(2)
  },
  downloadFab: {
    zIndex: 1,
    right: 32,
    bottom: 32,
    position: 'fixed !important'
  }
}))

export default function ArchiveEntryView(props) {
  const classes = useStyles()
  const {entryId} = props
  const {api} = useApi()
  const {raiseError} = useErrors()

  const [apiData, setApiData] = useState()
  const [doesNotExist, setDoesNotExist] = useState(false)

  useEffect(() => {
    api.get(`/entries/${entryId}/archive`, null, {returnRequest: true, jsonResponse: true})
      .then(apiData => {
        setApiData(apiData)
      })
      .catch(error => {
        if (error.name === 'DoesNotExist') {
          setDoesNotExist(true)
        } else {
          raiseError(error)
        }
      })
  }, [setApiData, setDoesNotExist, api, raiseError, entryId])

  const data = useMemo(() => {
    const archive = apiData?.response?.data?.archive
    if (!archive) {
      return null
    }
    const cleanedArchive = {...archive}
    cleanedArchive.processing_logs = undefined
    return cleanedArchive
  }, [apiData])

  if (doesNotExist) {
    return (
      <Page>
        <Typography className={classes.error}>
          No archive exists for this entry. Either the archive was not generated due
          to parsing or other processing errors (check the log tab), or the entry it
          self does not exist.
        </Typography>
      </Page>
    )
  }

  return (
    <Page width={'100%'}>
      {
        data && typeof data !== 'string'
          ? <div className={classes.archiveBrowser}>
            <ApiDataContext.Provider value={apiData}>
              <ArchiveBrowser data={data} />
            </ApiDataContext.Provider>
          </div> : <div>{
            data
              ? <div>
                <Typography>Processed data is not valid JSON. Displaying plain text instead.</Typography>
                <Card>
                  <CardContent>
                    <pre>{data || ''}</pre>
                  </CardContent>
                </Card>
              </div>
              : <Typography>loading ...</Typography>
          }</div>
      }
    </Page>
  )
}
ArchiveEntryView.propTypes = {
  entryId: PropTypes.string.isRequired
}
