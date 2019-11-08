import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { withApi } from './api'
import Search from './search/Search'
import SearchContext from './search/SearchContext'
import { Typography, Link } from '@material-ui/core'

export const help = `
This page allows you to **inspect** and **download** NOMAD datasets. It alsow allows you
to explore a dataset with similar controls that the search page offers.
`

class DatasetPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    datasetId: PropTypes.string.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    description: {
      padding: theme.spacing.unit * 3
    }
  })

  state = {
    dataset: {}
  }

  update() {
    const {datasetId, raiseError, api} = this.props
    api.search({
      owner: 'all',
      dataset_id: datasetId,
      page: 1, per_page: 1
    }).then(data => {
      const entry = data.results[0]
      const dataset = entry ? entry.datasets.find(ds => ds.id + '' === datasetId) : {}
      this.setState({dataset: dataset || {}})
    }).catch(error => {
        this.setState({dataset: {}})
        raiseError(error)
    })
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api || prevProps.datasetId !== this.props.datasetId) {
      this.update()
    }
  }

  render() {
    const { classes, datasetId } = this.props
    const { dataset } = this.state

    return (
      <div>
        <div className={classes.description}>
          <Typography variant="h4">{dataset.name || 'loading ...'}</Typography>
          <Typography>
            dataset{dataset.doi ? <span>, with DOI <Link href={dataset.doi}>{dataset.doi}</Link></span> : ''}
          </Typography>
        </div>

        <SearchContext query={{dataset_id: datasetId}} ownerTypes={['all', 'public']} >
          <Search resultTab="entries"/>
        </SearchContext>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(DatasetPage.styles))(DatasetPage)
