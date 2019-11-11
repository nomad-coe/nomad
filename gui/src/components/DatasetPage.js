import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { withApi, DoesNotExist } from './api'
import Search from './search/Search'
import SearchContext from './search/SearchContext'
import { Typography, Link } from '@material-ui/core'
import { DatasetActions } from './search/DatasetList'
import { withRouter } from 'react-router'

export const help = `
This page allows you to **inspect** and **download** NOMAD datasets. It alsow allows you
to explore a dataset with similar controls that the search page offers.
`

class DatasetPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    datasetId: PropTypes.string.isRequired,
    raiseError: PropTypes.func.isRequired,
    history: PropTypes.object.isRequired
  }

  static styles = theme => ({
    description: {
      flexGrow: 1,
      marginRight: theme.spacing.unit
    },
    header: {
      display: 'flex',
      flexDirection: 'row',
      padding: theme.spacing.unit * 3,
    },
    actions: {}
  })

  constructor(props) {
    super(props)
    this.handleChange = this.handleChange.bind(this)
  }

  state = {
    dataset: {},
    update: 0
  }

  update() {
    const {datasetId, raiseError, api} = this.props
    api.search({
      owner: 'all',
      dataset_id: datasetId,
      page: 1, per_page: 1
    }).then(data => {
      const entry = data.results[0]
      const dataset = entry && entry.datasets.find(ds => ds.id + '' === datasetId)
      if (!dataset) {
        this.setState({dataset: {}})
        raiseError(new DoesNotExist('Dataset does not exist any more or is not visible to you.'))
      }
      this.setState({dataset: {
        ...dataset, example: entry
      }})
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
      this.setState({dataset: {}}, () => this.update())
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
    const { classes, datasetId } = this.props
    const { dataset, update } = this.state

    return (
      <div>
        <div className={classes.header}>
          <div className={classes.description}>
            <Typography variant="h4">{dataset.name || 'loading ...'}</Typography>
            <Typography>
              dataset{dataset.doi ? <span>, with DOI <Link href={dataset.doi}>{dataset.doi}</Link></span> : ''}
            </Typography>
          </div>

          <div className={classes.actions}>
            {dataset && dataset.example && <DatasetActions
              dataset={dataset}
              onChange={this.handleChange}/>
            }
          </div>
        </div>

        <SearchContext query={{dataset_id: datasetId}} ownerTypes={['all', 'public']} update={update}>
          <Search resultTab="entries"/>
        </SearchContext>
      </div>
    )
  }
}

export default compose(withRouter, withApi(false), withErrors, withStyles(DatasetPage.styles))(DatasetPage)
