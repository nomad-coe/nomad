import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { withApi } from './api'
import Search from './search/Search'
import { Typography, Link, Fab } from '@material-ui/core'
import Download from './entry/Download'
import DownloadIcon from '@material-ui/icons/CloudDownload'

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
    root: {
    },
    description: {
      padding: theme.spacing.unit * 3
    },
    downloadFab: {
      zIndex: 1,
      right: 32,
      top: 56 + 32,
      position: 'fixed !important'
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
      <div className={classes.root}>
        <div className={classes.description}>
          <Typography variant="h4">{dataset.name || 'loading ...'}</Typography>
          <Typography>
            dataset{dataset.doi ? <span>, with DOI <Link href={dataset.doi}>{dataset.doi}</Link></span> : ''}
          </Typography>
        </div>
        <Search searchParameters={{owner: 'all', dataset_id: datasetId}} />
        <Download
          classes={{root: classes.downloadFab}} tooltip="download all rawfiles"
          component={Fab} className={classes.downloadFab} color="primary" size="medium"
          url={`raw/query?dataset_id=${datasetId}`} fileName={`${dataset.name}.json`}
        >
          <DownloadIcon />
        </Download>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(DatasetPage.styles))(DatasetPage)
