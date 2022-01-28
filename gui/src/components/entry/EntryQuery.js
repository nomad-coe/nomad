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
import { Typography, Link, makeStyles } from '@material-ui/core'
import { useHistory, useLocation } from 'react-router'
import qs from 'qs'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { getUrl } from '../nav/Routes'

const useStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(3)
  }
}))

export default function EntryQuery(props) {
  const classes = useStyles()
  const history = useHistory()
  const location = useLocation()
  const {api} = useApi()
  const {raiseError} = useErrors()

  const [doesNotExist, setDoesNotExist] = useState(false)

  const queryParams = useMemo(() => qs.parse(location.search.substring(1)), [location])

  useEffect(() => {
    api.get('/entries', {...queryParams})
      .then(response => {
        if (response.pagination.total > 0) {
          const {entry_id, upload_id} = response.data[0]
          history.push(getUrl(`entry/id/${upload_id}/${entry_id}`))
        } else {
          setDoesNotExist(true)
        }
      })
      .catch(raiseError)
  }, [history, api, raiseError, setDoesNotExist, queryParams])

  let message = 'loading ...'

  if (doesNotExist) {
    if (queryParams && queryParams['external_id'] && queryParams['external_id'].startsWith('mp-')) {
      message = <React.Fragment>
        This particular entry <Link href={`https://materialsproject.org/tasks/${queryParams['external_id']}#`}>
          {queryParams['external_id']}
        </Link> has not yet been provided to NOMAD by the Materials Project.
      </React.Fragment>
    } else if (api.isLoggedIn) {
      message = `
          This URL points to an entry that either does not exist, or that you are not
          authorized to see.`
    } else {
      message = `
          This URL points to an entry that either does not exist, or that is not
          publically visibile. Please login; you might be authorized to view it.`
    }
  }

  return (
    <Typography className={classes.root}>{message}</Typography>
  )
}
