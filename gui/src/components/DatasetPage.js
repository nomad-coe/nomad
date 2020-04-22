import React, { useContext, useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withErrors, errorContext } from './errors'
import { withApi, apiContext } from './api'
import Search from './search/Search'
import { Typography, makeStyles } from '@material-ui/core'
import { DatasetActions, DOI } from './search/DatasetList'
import { matchPath, useLocation, useHistory, useRouteMatch } from 'react-router'

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
  const history = useHistory()

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
  }, [location.pathname, api])

  const handleChange = dataset => {
    if (dataset) {
      setDataset({dataset: dataset})
    } else {
      history.goBack()
    }
  }

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

      <div className={classes.actions}>
        {dataset && dataset.example && <DatasetActions
          dataset={dataset}
          onChange={handleChange}/>
        }
      </div>
    </div>

    <Search
      initialQuery={{owner: 'all'}}
      query={{dataset_id: datasetId}}
      ownerTypes={['all', 'public']}
      initialResultTab="entries" availableResultTabs={['entries', 'groups', 'datasets']}
    />
  </div>
}

DatasetPage.propTypes = {
  classes: PropTypes.object.isRequired,
  api: PropTypes.object.isRequired,
  raiseError: PropTypes.func.isRequired,
  location: PropTypes.object.isRequired,
  match: PropTypes.object.isRequired,
  history: PropTypes.object.isRequired
}
