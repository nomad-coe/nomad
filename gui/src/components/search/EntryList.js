import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Link, Typography, Tooltip, IconButton, TablePagination } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import { withDomain } from '../domains'
import DataTable from '../DataTable'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import DetailsIcon from '@material-ui/icons/MoreHoriz'
import EditUserMetadataDialog from '../EditUserMetadataDialog'
import DownloadButton from '../DownloadButton'

export class EntryListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    query: PropTypes.object.isRequired,
    onChange: PropTypes.func,
    history: PropTypes.any.isRequired,
    order_by: PropTypes.string.isRequired,
    order: PropTypes.number.isRequired,
    page: PropTypes.number.isRequired,
    per_page: PropTypes.number.isRequired,
    domain: PropTypes.object.isRequired,
    editable: PropTypes.bool,
    columns: PropTypes.object,
    title: PropTypes.string,
    actions: PropTypes.element,
    selectedColumns: PropTypes.arrayOf(PropTypes.string)
  }

  static styles = theme => ({
    root: {
      overflow: 'auto'
    },
    entryDetails: {
      display: 'flex'
    },
    entryDetailsRow: {
      paddingRight: theme.spacing.unit * 3
    }
  })

  state = {
    selected: []
  }

  static defaultColumns = {
    mainfile: {
      label: 'Mainfile',
      render: entry => entry.mainfile,
      supportsSort: true,
      description: 'The mainfile of this entry.'
    },
    upload_time: {
      label: 'Upload time',
      render: entry => new Date(entry.upload_time).toLocaleString(),
      supportsSort: true,
      description: 'The time this entry was uploaded.'
    },
    authors: {
      label: 'Authors',
      render: entry => entry.authors.map(author => author.name).join('; '),
      supportsSort: true,
      description: 'The authors of this entry. This includes the uploader and its co-authors.'
    },
    co_authors: {
      label: 'co-Authors',
      render: entry => entry.authors.filter(user => user.user_id !== entry.uploader.user_id).map(author => author.name).join('; '),
      supportsSort: false,
      description: 'The people that this entry was co authored with'
    },
    shared_with: {
      label: 'Shared with',
      render: entry => entry.owners.filter(user => user.user_id !== entry.uploader.user_id).map(author => author.name).join('; '),
      supportsSort: false,
      description: 'The people that this entry was shared with'
    },
    uploader: {
      label: 'Uploader',
      render: entry => entry.uploader.name,
      supportsSort: true,
      description: 'The uploader of this entry.'
    },
    comment: {
      label: 'Comment',
      render: entry => entry.comment,
      supportsSort: false,
      description: 'User provided comment on this entry'
    },
    references: {
      label: 'References',
      render: entry => {
        const refs = entry.references || []
        if (refs.length > 0) {
          return (
            <div style={{display: 'inline'}}>
              {refs.map((ref, i) => <span key={ref}>
                <a href={ref}>{ref}</a>{(i + 1) < refs.length ? ', ' : <React.Fragment/>}
              </span>)}
            </div>
          )
        } else {
          return <i>no references</i>
        }
      },
      supportsSort: true
    },
    datasets: {
      label: 'Datasets',
      render: entry => {
        const datasets = entry.datasets || []
        if (datasets.length > 0) {
          return datasets.map(dataset => dataset.name).join(', ')
        } else {
          return <i>no datasets</i>
        }
      },
      supportsSort: false,
      description: 'The dataset names that this entry belongs to.'
    }
  }

  componentWillUpdate(prevProps) {
    if (prevProps.data !== this.props.data) {
      this.setState({selected: []})
    }
  }

  handleChange(changes) {
    if (this.props.onChange) {
      this.props.onChange(changes)
    }
  }

  handleClickCalc(calc) {
    this.props.history.push(`/entry/id/${calc.upload_id}/${calc.calc_id}`)
  }

  handleChangePage = (event, page) => {
    this.handleChange({page: page + 1})
  }

  handleChangeRowsPerPage = event => {
    const rowsPerPage = event.target.value
    this.handleChange({per_page: rowsPerPage})
  }

  handleSort(columnKey) {
    if (this.props.order_by === columnKey) {
      this.handleChange({order: this.props.order * -1})
    } else {
      this.handleChange({order: 1, order_by: columnKey})
    }
  }

  selectionQuery() {
    const { selected } = this.state
    if (selected) {
      return {
        'calc_id': selected.join(',')
      }
    } else {
      return this.props.query
    }
  }

  renderEntryDetails(row) {
    const { classes, domain } = this.props
    return (<div className={classes.entryDetails}>
      <div className={classes.entryDetailsRow}>
        <Quantity column>
          <Quantity row>
            <Quantity quantity="formula" label='formula' noWrap data={row} />
          </Quantity>
          <Quantity row>
            <Quantity quantity="code_name" label='dft code' noWrap data={row} />
            <Quantity quantity="code_version" label='dft code version' noWrap data={row} />
          </Quantity>
          <Quantity row>
            <Quantity quantity="basis_set" label='basis set' noWrap data={row} />
            <Quantity quantity="xc_functional" label='xc functional' noWrap data={row} />
          </Quantity>
          <Quantity row>
            <Quantity quantity="system" label='system type' noWrap data={row} />
            <Quantity quantity="crystal_system" label='crystal system' noWrap data={row} />
            <Quantity quantity='spacegroup_symbol' label="spacegroup" noWrap data={row} >
              <Typography noWrap>
                {row.spacegroup_symbol} ({row.spacegroup})
              </Typography>
            </Quantity>
          </Quantity>
        </Quantity>
      </div>

      <div className={classes.entryDetailsRow} style={{flexGrow: 1}}>
        <Quantity className={classes.entryDetailsRow} column>
          <Quantity quantity='comment' placeholder='no comment' data={row} />
          <Quantity quantity='references' placeholder='no references' data={row}>
            <div>
              {(row.references || []).map(ref => <Typography key={ref} noWrap>
                <a href={ref}>{ref}</a>
              </Typography>)}
            </div>
          </Quantity>
          <Quantity quantity='authors' data={row}>
            <Typography>
              {(row.authors || []).map(author => author.name).join('; ')}
            </Typography>
          </Quantity>
          <Quantity quantity='datasets' placeholder='no datasets' data={row}>
            <div>
              {(row.datasets || []).map(ds => (
                <Typography key={ds.id}>
                  <Link component={RouterLink} to={`/dataset/id/${ds.id}`}>{ds.name}</Link>
                  {ds.doi ? <span>&nbsp; (<Link href={ds.doi}>{ds.doi}</Link>)</span> : <React.Fragment/>}
                </Typography>))}
            </div>
          </Quantity>
        </Quantity>
      </div>

      <div className={classes.entryDetailsRow} style={{maxWidth: '33%'}}>
        <Quantity column >
          {/* <Quantity quantity="pid" label='PID' placeholder="not yet assigned" noWrap data={row} withClipboard /> */}
          <Quantity quantity="upload_id" label='upload id' data={row} noWrap withClipboard />
          <Quantity quantity="calc_id" label={`${domain.entryLabel} id`} noWrap withClipboard data={row} />
          <Quantity quantity='mainfile' noWrap data={row} withClipboard />
          <Quantity quantity="upload_time" label='upload time' noWrap data={row} >
            <Typography noWrap>
              {new Date(row.upload_time * 1000).toLocaleString()}
            </Typography>
          </Quantity>
          {/* <Quantity quantity="calc_hash" label={`${domain.entryLabel} hash`} noWrap data={row} />
          <Quantity quantity="raw_id" label='raw id' noWrap data={row} withClipboard />
          <Quantity quantity="last_processing" label='last processing' placeholder="not processed" noWrap data={row}>
            <Typography noWrap>
              {new Date(row.last_processing * 1000).toLocaleString()}
            </Typography>
          </Quantity> */}
          <Quantity quantity="last_processing" label='processing version' noWrap placeholder="not processed" data={row}>
            <Typography noWrap>
              {row.nomad_version}/{row.nomad_commit}
            </Typography>
          </Quantity>
        </Quantity>
      </div>
    </div>)
  }

  handleViewEntryPage(event, row) {
    event.stopPropagation()
    this.props.history.push(`/entry/id/${row.upload_id}/${row.calc_id}`)
  }

  renderEntryActions(row, selected) {
    return <Tooltip title="View entry page">
      <IconButton style={selected ? {color: 'white'} : null} onClick={event => this.handleViewEntryPage(event, row)}>
        <DetailsIcon />
      </IconButton>
    </Tooltip>
  }

  render() {
    const { classes, data, order, order_by, page, per_page, domain, editable, title, query, actions, ...rest } = this.props
    const { selected } = this.state

    const results = data.results || []
    const total = data.pagination && data.pagination.total
    const totalNumber = total || 0

    const columns = this.props.columns || {
      ...domain.searchResultColumns,
      ...EntryListUnstyled.defaultColumns
    }

    const defaultSelectedColumns = this.props.selectedColumns || [
      ...domain.defaultSearchResultColumns,
      'datasets', 'authors']

    const pagination = <TablePagination
      count={totalNumber}
      rowsPerPage={per_page}
      page={page - 1}
      onChangePage={this.handleChangePage}
      onChangeRowsPerPage={this.handleChangeRowsPerPage}
      labelDisplayedRows={({ from, to, count }) => `${from.toLocaleString()}-${to.toLocaleString()} of ${count.toLocaleString()}`}
    />

    const example = selected && selected.length > 0 ? results.find(d => d.calc_id === selected[0]) : results[0]
    const selectQuery = selected ? {calc_id: selected.join(',')} : query
    const createActions = (props, moreActions) => <React.Fragment>
      {example && editable ? <EditUserMetadataDialog
        example={example} total={totalNumber}
        onEditComplete={() => this.props.onChange()}
        {...props}
      /> : ''}
      <DownloadButton
        tooltip="Download raw files"
        {...props}/>
      {moreActions}
    </React.Fragment>
    const selectActions = createActions({query: selectQuery, buttonProps: {color: 'secondary'}})
    const allActions = actions

    return (
      <div className={classes.root}>
        <DataTable
          entityLabels={['domain.entryLabel', domain.entryLabel + 's']}
          selectActions={selectActions}
          id={row => row.calc_id}
          total={total}
          columns={columns}
          selectedColumns={defaultSelectedColumns}
          entryDetails={this.renderEntryDetails.bind(this)}
          entryActions={this.renderEntryActions.bind(this)}
          data={results}
          order={order === 1 ? 'desc' : 'asc'}
          orderBy={order_by}
          selected={this.state.selected}
          onSelectionChanged={selection => this.setState({selected: selection})}
          onOrderChanged={(order, orderBy) => this.handleChange({order: order === 'asc' ? -1 : 1, order_by: orderBy})}
          rows={per_page}
          pagination={pagination}
          actions={allActions}
          {...rest}
        />
      </div>
    )
  }
}

const EntryList = compose(withRouter, withDomain, withStyles(EntryListUnstyled.styles))(EntryListUnstyled)

export default EntryList
