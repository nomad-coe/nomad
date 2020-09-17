import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, TableCell, Toolbar, IconButton, FormGroup, Tooltip, Link } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import NextIcon from '@material-ui/icons/ChevronRight'
import StartIcon from '@material-ui/icons/SkipPrevious'
import DataTable from '../DataTable'
import SearchIcon from '@material-ui/icons/Search'
import DOIIcon from '@material-ui/icons/Bookmark'
import DeleteIcon from '@material-ui/icons/Delete'
import { withApi } from '../api'
import EditUserMetadataDialog from '../EditUserMetadataDialog'
import DownloadButton from '../DownloadButton'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import ConfirmDialog from '../uploads/ConfirmDialog'
import { oasis } from '../../config'
import { authorList } from '../../utils'

class DOIUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    doi: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {
      display: 'inline-flex',
      alignItems: 'center',
      flexDirection: 'row',
      flexWrap: 'nowrap'
    }
  })

  render() {
    const {classes, doi} = this.props
    const url = `https://dx.doi.org/${doi}`
    return <span className={classes.root}>
      <Link href={url}>{doi}</Link>
      <CopyToClipboard
        text={url} onCopy={() => null}
      >
        <Tooltip title={`Copy DOI to clipboard`}>
          <IconButton style={{margin: 3, marginRight: 0, padding: 4}}>
            <ClipboardIcon style={{fontSize: 16}} />
          </IconButton>
        </Tooltip>
      </CopyToClipboard>
    </span>
  }
}

export const DOI = withStyles(DOIUnstyled.styles)(DOIUnstyled)

class DatasetActionsUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    dataset: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired,
    search: PropTypes.bool,
    user: PropTypes.object,
    onChange: PropTypes.func,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    group: {
      flexWrap: 'nowrap',
      flexDirection: 'row-reverse'
    }
  })

  state = {
    confirmDoi: false
  }

  constructor(props) {
    super(props)
    this.handleClickDOI = this.handleClickDOI.bind(this)
    this.handleClickDataset = this.handleClickDataset.bind(this)
    this.handleClickDelete = this.handleClickDelete.bind(this)
    this.handleEdit = this.handleEdit.bind(this)
  }

  handleClickDataset() {
    const {dataset: {id}} = this.props
    this.props.history.push(`/dataset/id/${id}`)
  }

  handleClickDOI(after) {
    const {api, dataset, onChange, raiseError} = this.props
    const datasetName = dataset.name

    api.assignDatasetDOI(datasetName)
      .then(dataset => {
        if (onChange) {
          onChange(dataset)
        }
        if (after) {
          after()
        }
      })
      .catch(raiseError)
  }

  handleEdit() {
    const {onChange, dataset} = this.props
    if (onChange) {
      onChange(dataset)
    }
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
      dataset.example.owners.find(author => author.user_id === user.sub)

    const canAssignDOI = !doi
    const canDelete = !doi
    const query = {dataset_id: [dataset.dataset_id]}

    return <FormGroup row classes={{root: classes.group}}>
      {search && <Tooltip title="Open a search page with entries from this dataset only.">
        <IconButton onClick={this.handleClickDataset}>
          <SearchIcon />
        </IconButton>
      </Tooltip>}
      {<DownloadButton query={query} tooltip="Download dataset" />}
      {editable && canDelete && <Tooltip title="Delete this dataset.">
        <IconButton onClick={this.handleClickDelete}>
          <DeleteIcon />
        </IconButton>
      </Tooltip>}
      {editable && <EditUserMetadataDialog
        title="Edit metadata of all dataset entries"
        example={dataset.example} query={query}
        total={dataset.total} onEditComplete={this.handleEdit}
      />}
      {!oasis && editable && canAssignDOI && <Tooltip title="Assign a DOI to this dataset.">
        <IconButton onClick={() => this.setState({confirmDoi: true})}>
          <DOIIcon />
        </IconButton>
      </Tooltip>}
      <ConfirmDialog
        open={this.state.confirmDoi}
        title="Assign a DOI"
        content={`
          DOIs are **permanent**. Are you sure that you want to assign a DOI to this
          dataset? Once the DOI was assigned, entries cannot removed from the dataset and
          the dataset cannot be deleted.
        `}
        onClose={() => this.setState({confirmDoi: false})}
        onConfirm={() => {
          this.handleClickDOI(() => this.setState({confirmDoi: false}))
        }}
      />
    </FormGroup>
  }
}

export const DatasetActions = compose(withRouter, withApi(false), withStyles(DatasetActionsUnstyled.styles))(DatasetActionsUnstyled)

class DatasetListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object,
    total: PropTypes.number,
    onChange: PropTypes.func.isRequired,
    onEdit: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired,
    datasets_after: PropTypes.string,
    per_page: PropTypes.number,
    actions: PropTypes.element
  }

  static styles = theme => ({
    root: {
      overflow: 'auto',
      paddingLeft: theme.spacing(2),
      paddingRight: theme.spacing(2)
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
      description: 'The name given to this dataset by its creator',
      render: (dataset) => dataset.name
    },
    created: {
      label: 'Created',
      description: 'The data when this dataset was created',
      render: (dataset) => dataset.created && new Date(dataset.created).toLocaleString()
    },
    DOI: {
      label: 'Dataset DOI',
      description: 'The DOI of the dataset, if a DOI was assigned',
      render: (dataset) => dataset.doi && <DOI doi={dataset.doi} />
    },
    entries: {
      label: 'Entries',
      description: 'Number of entries that comprise the group',
      render: (dataset) => dataset.total
    },
    authors: {
      label: 'Authors',
      description: 'Authors including the uploader and the co-authors',
      render: (dataset) => {
        const authors = dataset.example.authors
        if (authors.length > 3) {
          return authorList(authors.filter((_, index) => index < 2)) + ' et al'
        } else {
          return authorList(authors)
        }
      }
    }
  }

  renderEntryActions(entry) {
    const {onEdit} = this.props
    return <DatasetActions search dataset={entry} onChange={onEdit} />
  }

  render() {
    const { classes, data, total, datasets_after, per_page, onChange, actions } = this.props
    const datasets = data.datasets_grouped || {values: []}
    const results = Object.keys(datasets.values).map(id => {
      const exampleDataset = datasets.values[id].examples[0].datasets.find(ds => ds.dataset_id === id)
      return {
        ...exampleDataset,
        id: id,
        total: datasets.values[id].total,
        example: datasets.values[id].examples[0]
      }
    })
    const after = datasets.after
    const perPage = per_page || 10

    let paginationText
    if (datasets_after) {
      paginationText = `next ${results.length.toLocaleString()} of ${(total || 0).toLocaleString()}`
    } else {
      paginationText = `1-${results.length.toLocaleString()} of ${(total || 0).toLocaleString()}`
    }

    const pagination = <TableCell colSpan={1000} classes={{root: classes.scrollCell}}>
      <Toolbar className={classes.scrollBar}>
        <span className={classes.scrollSpacer}>&nbsp;</span>
        <span>{paginationText}</span>
        <IconButton disabled={!datasets_after} onClick={() => onChange({datasets_grouped_after: null})}>
          <StartIcon />
        </IconButton>
        <IconButton disabled={results.length < perPage} onClick={() => onChange({datasets_grouped_after: after})}>
          <NextIcon />
        </IconButton>
      </Toolbar>
    </TableCell>

    return <DataTable
      entityLabels={['dataset', 'datasets']}
      id={row => row.id}
      total={total}
      columns={this.columns}
      selectedColumns={['name', 'DOI', 'entries', 'authors']}
      selectedColumnsKey="datasets"
      entryActions={this.renderEntryActions}
      data={results}
      rows={perPage}
      actions={actions}
      pagination={pagination}
    />
  }
}

const DatasetList = compose(withRouter, withApi(false), withStyles(DatasetListUnstyled.styles))(DatasetListUnstyled)

export default DatasetList
