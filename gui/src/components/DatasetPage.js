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
import { errorContext } from './errors'
import { apiContext } from './api'
import Search from './search/Search'
import { Typography, makeStyles } from '@material-ui/core'
import { useLocation, useRouteMatch } from 'react-router'
import { DOI } from './search/results/DatasetList'

export const help = `
This page allows you to **inspect** and **download** NOMAD datasets. It also allows you
to explore a dataset with similar controls that the search page offers.
`

const useStyles = makeStyles(theme => ({
  description: {
    flexGrow: 1,
    marginRight: theme.spacing(1)
  },
  header: {
    display: 'flex',
    flexDirection: 'row',
    padding: theme.spacing(3)
  }
}))

export default function DatasetPage() {
  const classes = useStyles()
  const [dataset, setDataset] = useState({})

  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const location = useLocation()
  const match = useRouteMatch()

  const {datasetId} = match.params

  useEffect(() => {
    api.search({
      owner: 'all',
      dataset_id: datasetId,
      page: 1,
      per_page: 1
    }).then(data => {
      const entry = data.results[0]
      const dataset = entry && entry.datasets.find(ds => ds.dataset_id + '' === datasetId)
      if (!dataset) {
        setDataset({isEmpty: true})
      }
      setDataset({...dataset, example: entry})
    }).catch(error => {
      setDataset({})
      raiseError(error)
    })
  }, [datasetId, location.pathname, api, raiseError])

  if (!dataset) {
    return <div>loading...</div>
  }

  return <div>
    <div className={classes.header}>
      <div className={classes.description}>
        <Typography variant="h4">{dataset.name || (dataset.isEmpty && 'Empty or non existing dataset') || 'loading ...'}</Typography>
        <Typography>
          dataset{dataset.doi ? <span>, with DOI <DOI doi={dataset.doi} /></span> : ''}
        </Typography>
      </div>
    </div>

    <Search
      initialQuery={{owner: 'all'}}
      query={{dataset_id: [datasetId]}}
      ownerTypes={['all', 'public']}
      initialResultTab="entries"
      resultListProps={{showAccessColumn: true}}
      availableResultTabs={['entries', 'groups', 'datasets']}
    />
  </div>
}
