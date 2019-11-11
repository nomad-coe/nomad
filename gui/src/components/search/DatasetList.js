import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, TableCell, Toolbar, IconButton, Button, FormGroup, Tooltip } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import { withDomain } from '../domains'
import NextIcon from '@material-ui/icons/ChevronRight'
import StartIcon from '@material-ui/icons/SkipPrevious'
import DataTable from '../DataTable'
import SearchIcon from '@material-ui/icons/Search'
import DOIIcon from '@material-ui/icons/Bookmark'
import DeleteIcon from '@material-ui/icons/Delete'
import { withApi } from '../api'

class DatasetActionsUnstyled extends React.Component {

  static propTypes = {
    classes: PropTypes.object.isRequired,
    dataset: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired,
    search: PropTypes.bool,
    user: PropTypes.object,
    onChange: PropTypes.func,
    api: PropTypes.object.isRequired
  }

  static styles = theme => ({
    group: {
      flexWrap: 'nowrap',
      flexDirection: 'row-reverse'
    }
  })

  constructor(props) {
    super(props)
    this.handleClickDOI = this.handleClickDOI.bind(this)
    this.handleClickDataset = this.handleClickDataset.bind(this)
    this.handleClickDelete = this.handleClickDelete.bind(this)
  }

  handleClickDataset() {
    const {dataset: {id}} = this.props
    this.props.history.push(`/dataset/id/${id}`)
  }

  handleClickDOI() {
    const {api, dataset, onChange, raiseError} = this.props
    const datasetName = dataset.name

    api.assignDatasetDOI(datasetName)
      .then(dataset => {
        if (onChange) {
          onChange(dataset)
        }
      })
      .catch(raiseError)
  }

  handleClickDelete() {
    const {api, dataset, onChange, raiseError} = this.props
    const datasetName = dataset.name

    api.deleteDataset(datasetName)
      .then(dataset => {
        if (onChange) {
          onChange(null)
        }
      })
      .catch(raiseError)
  }


  render() {
    const {dataset, search, user, classes} = this.props
    const {doi} = dataset
    const editable = user && dataset.example &&
      dataset.example.authors.find(author => author.user_id === user.sub)

    const canAssignDOI = !doi
    const canDelete = !doi

    return <FormGroup row classes={{root: classes.group}}>
      {search && <Tooltip title="Open a search page with entries from this dataset only.">
        <IconButton onClick={this.handleClickDataset}>
          <SearchIcon />
        </IconButton>
      </Tooltip>}
      {editable && canDelete && <Tooltip title="Delete this dataset.">
        <IconButton onClick={this.handleClickDelete}>
          <DeleteIcon />
        </IconButton>
      </Tooltip>}
      {editable && canAssignDOI && <Tooltip title="Assign a DOI to this dataset.">
        <IconButton onClick={this.handleClickDOI}>
          <DOIIcon />
        </IconButton>
      </Tooltip>}
    </FormGroup>
  }
}

export const DatasetActions = compose(withRouter, withApi(false), withStyles(DatasetActionsUnstyled.styles))(DatasetActionsUnstyled)

class DatasetListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    total: PropTypes.number.isRequired,
    onChange: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired,
    datasets_after: PropTypes.string,
    actions: PropTypes.element
  }

  static styles = theme => ({
    root: {
      overflow: 'auto',
      paddingLeft: theme.spacing.unit * 2,
      paddingRight: theme.spacing.unit * 2
    },
    scrollCell: {
      padding: 0
    },
    scrollBar: {
      minHeight: 56,
      padding: 0
    },
    scrollSpacer: {
      flexGrow: 1
    },
    clickableRow: {
      cursor: 'pointer'
    }
  })

  constructor(props) {
    super(props)
    this.renderEntryActions = this.renderEntryActions.bind(this)
  }

  columns = {
    name: {
      label: 'Dataset name',
      render: (dataset) => dataset.name
    },
    DOI: {
      label: 'Dataset DOI',
      render: (dataset) => dataset.doi
    },
    entries: {
      label: 'Entries',
      render: (dataset) => dataset.total
    },
    authors: {
      label: 'Authors',
      render: (dataset) => {
        const authors = dataset.example.authors
        if (authors.length > 3) {
          return authors.filter((_, index) => index < 2).map(author => author.name).join('; ') + ' et al'
        } else {
          return authors.map(author => author.name).join('; ')
        }
      }
    }
  }

  renderEntryActions(entry) {
    const {onChange} = this.props
    return <DatasetActions search dataset={entry} onChange={() => onChange({})} />
  }

  render() {
    const { classes, data, total, datasets_after, onChange, actions } = this.props
    const datasets = data.datasets || {values: []}
    const results = Object.keys(datasets.values).map(id => {
      const exampleDataset = datasets.values[id].examples[0].datasets.find(ds => ds.id === id)
      return {
        ...exampleDataset,
        id: id,
        total: datasets.values[id].total,
        example: datasets.values[id].examples[0]
      }
    })
    const per_page = 10
    const after = datasets.after

    let paginationText
    if (datasets_after) {
      paginationText = `next ${results.length} of ${total}`
    } else {
      paginationText = `1-${results.length} of ${total}`
    }

    const pagination = <TableCell colSpan={1000} classes={{root: classes.scrollCell}}>
      <Toolbar className={classes.scrollBar}>
        <span className={classes.scrollSpacer}>&nbsp;</span>
        <span>{paginationText}</span>
        <IconButton disabled={!datasets_after} onClick={() => onChange({datasets_after: null})}>
          <StartIcon />
        </IconButton>
        <IconButton disabled={results.length < per_page} onClick={() => onChange({datasets_after: after})}>
          <NextIcon />
        </IconButton>
      </Toolbar>
    </TableCell>

    return <DataTable
      title={`${total.toLocaleString()} datasets`}
      id={row => row.id}
      total={total}
      columns={this.columns}
      // selectedColumns={defaultSelectedColumns}
      // entryDetails={this.renderEntryDetails.bind(this)}
      entryActions={this.renderEntryActions}
      data={results}
      rows={per_page}
      actions={actions}
      pagination={pagination}
    />
  }
}

const DatasetList = compose(withRouter, withDomain, withApi(false), withStyles(DatasetListUnstyled.styles))(DatasetListUnstyled)

export default DatasetList
