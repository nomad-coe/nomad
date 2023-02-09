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
import React, { useState, useEffect } from 'react'
import { Typography, makeStyles } from '@material-ui/core'
import { useErrors } from '../errors'
import FileBrowser from '../archive/FileBrowser'
import { useApi } from '../api'
import Page from '../Page'
import { useEntryStore } from './EntryContext'
import { createUploadUrl } from '../../utils'
import { apiBase } from '../../config'

const useStyles = makeStyles(theme => ({
  error: {
    marginTop: theme.spacing(2)
  }
}))
const BrowseEntryFilesView = React.memo((props) => {
  const {entryId} = useEntryStore()
  const classes = useStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()

  const [data, setData] = useState(null)
  const [doesNotExist, setDoesNotExist] = useState(false)

  useEffect(() => {
    api.get(`/entries/${entryId}?include=mainfile&include=upload_id&include=writers&include=published`)
      .then(response => {
        setData(response.data)
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
        <Typography className={classes.error}>
          This entry does not exist.
        </Typography>
      </Page>
    )
  }
  if (data && data?.upload_id && data?.mainfile) {
    const lastSlash = data.mainfile.lastIndexOf('/')
    const mainfileDirname = lastSlash === -1 ? '' : data.mainfile.substr(0, lastSlash)
    const mainfileBasename = data.mainfile.split('/').pop()
    return <Page>
      <FileBrowser
        uploadUrl={createUploadUrl(apiBase, data.upload_id, mainfileDirname)} // TODO: deploymentUrl should come from entry context
        rootTitle="Entry files"
        highlightedItem={mainfileBasename}
      />
    </Page>
  } else {
    return <Page><Typography>loading ...</Typography></Page>
  }
})
BrowseEntryFilesView.propTypes = {}
export default BrowseEntryFilesView
