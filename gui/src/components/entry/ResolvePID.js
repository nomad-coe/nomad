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
import { Typography, makeStyles } from '@material-ui/core'
import { matchPath, useLocation, useRouteMatch, useHistory } from 'react-router'
import {useApi} from '../api'
import {useErrors} from '../errors'
import { getUrl } from '../nav/Routes'

const useStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(3)
  }
}))

export default function ResolvePID() {
  const classes = useStyles()
  const {api, user} = useApi()
  const {raiseError} = useErrors()
  const history = useHistory()
  const location = useLocation()
  const match = useRouteMatch()

  const [doesNotExist, setDoesNotExist] = useState(false)

  useEffect(() => {
    const pidMatch = matchPath(location.pathname, {
      path: `${match.path}/:pid/:handle?`
    })
    const { pid, handle } = pidMatch.params
    const pidWithHandle = handle ? pid + '/' + handle : pid

    api.post('/entries/query', {owner: 'all', query: {pid: pidWithHandle}})
      .then(response => {
        if (response.pagination.total >= 1) {
          const entry = response.data[0]
          history.push(getUrl(`entry/id/${entry.upload_id}/${entry.entry_id}`, location))
        } else {
          setDoesNotExist(true)
        }
      })
      .catch(raiseError)
  }, [setDoesNotExist, history, location, match, api, raiseError])

  let message = 'loading ...'

  if (doesNotExist) {
    if (user) {
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
