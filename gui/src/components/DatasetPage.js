import React, { useContext, useState, useEffect } from 'react'
import { errorContext } from './errors'
import { apiContext } from './api'
import Search from './search/Search'
import { Typography, makeStyles } from '@material-ui/core'
import { matchPath, useLocation, useRouteMatch } from 'react-router'
import { DOI } from './search/DatasetList'

export const help = `
This page allows you to **inspect** and **download** NOMAD datasets. It alsow allows you
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

  const {datasetId} = matchPath(location.pathname, {
    path: `${match.path}/:datasetId`
  }).params

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

  console.log('### DatasetPage', dataset)
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
      availableResultTabs={['entries', 'groups', 'datasets']}
    />
  </div>
}
