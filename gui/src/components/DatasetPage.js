import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { withApi } from './api'
import Search from './search/Search'
import SearchContext from './search/SearchContext'
import { Typography } from '@material-ui/core'
import { DatasetActions, DOI } from './search/DatasetList'
import { matchPath } from 'react-router'

export const help = `
This page allows you to **inspect** and **download** NOMAD datasets. It alsow allows you
to explore a dataset with similar controls that the search page offers.
`

class DatasetPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    location: PropTypes.object.isRequired,
    match: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired
  }

  static styles = theme => ({
    description: {
      flexGrow: 1,
      marginRight: theme.spacing(1)
    },
    header: {
      display: 'flex',
      flexDirection: 'row',
      padding: theme.spacing(3)
    },
    actions: {}
  })

  constructor(props) {
    super(props)
    this.handleChange = this.handleChange.bind(this)
  }

  state = {
    dataset: {},
    empty: false,
    update: 0
  }

  datasetId() {
    const { location, match } = this.props
    const pidMatch = matchPath(location.pathname, {
      path: `${match.path}/:datasetId`
    })
    let { datasetId } = pidMatch.params
    return datasetId
  }

  update() {
    const { api, raiseError } = this.props
    const datasetId = this.datasetId()
    api.search({
      owner: 'all',
      dataset_id: datasetId,
      page: 1,
      per_page: 1
    }).then(data => {
      const entry = data.results[0]
      const dataset = entry && entry.datasets.find(ds => ds.dataset_id + '' === datasetId)
      if (!dataset) {
        this.setState({dataset: {}, empty: true})
      }
      this.setState({dataset: {
        ...dataset, example: entry, empty: false
      }})
    }).catch(error => {
      this.setState({dataset: {}, empty: false})
      raiseError(error)
    })
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.location.pathname !== this.props.location.pathname || prevProps.api !== this.props.api) {
      this.setState({dataset: {}, empty: false}, () => this.update())
    }
  }

  handleChange(dataset) {
    if (dataset) {
      this.setState({dataset: dataset, update: this.state.update + 1})
    } else {
      this.props.history.goBack()
    }
  }

  render() {
    const { classes } = this.props
    const { dataset, update, empty } = this.state
    const datasetId = this.datasetId()

    return (
      <div>
        <div className={classes.header}>
          <div className={classes.description}>
            <Typography variant="h4">{dataset.name || (empty && 'Empty or non existing dataset') || 'loading ...'}</Typography>
            <Typography>
              dataset{dataset.doi ? <span>, with DOI <DOI doi={dataset.doi} /></span> : ''}
            </Typography>
          </div>

          <div className={classes.actions}>
            {dataset && dataset.example && <DatasetActions
              dataset={dataset}
              onChange={this.handleChange}/>
            }
          </div>
        </div>

        <SearchContext
          initialQuery={{owner: 'all'}}
          query={{dataset_id: datasetId}}
          ownerTypes={['all', 'public']} update={update}
        >
          <Search
            resultTab="entries" tabs={['entries', 'groups', 'datasets']}
          />
        </SearchContext>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(DatasetPage.styles))(DatasetPage)
