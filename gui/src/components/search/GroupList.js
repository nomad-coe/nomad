import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, TableCell, Toolbar, IconButton, Table, TableHead, TableRow, TableBody, Tooltip } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import NextIcon from '@material-ui/icons/ChevronRight'
import StartIcon from '@material-ui/icons/SkipPrevious'
import DataTable from '../DataTable'
import { withApi } from '../api'
import { EntryListUnstyled } from './EntryList'
import MoreIcon from '@material-ui/icons/MoreHoriz'
import DownloadButton from '../DownloadButton'
import SearchContext from './SearchContext'

class GroupUnstyled extends React.Component {
  static contextType = SearchContext.type

  static propTypes = {
    classes: PropTypes.object.isRequired,
    groupHash: PropTypes.string.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    history: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    }
  })

  state = {
    entries: []
  }

  update() {
    const {groupHash, api, raiseError} = this.props
    const {query} = this.context.state
    api.search({...query, 'dft.group_hash': groupHash, per_page: 100})
      .then(data => {
        this.setState({entries: data.results})
      })
      .catch(raiseError)
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.groupHash !== this.props.groupHash || prevProps.api !== this.props.api) {
      this.update()
    }
  }

  render() {
    const {history} = this.props
    const {entries} = this.state
    return (
      <Table>
        <TableHead>
          <TableRow>
            <TableCell>Mainfile</TableCell>
            <TableCell>Upload time</TableCell>
            <TableCell></TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {entries.map(entry => (
            <TableRow key={entry.calc_id}>
              <TableCell>{entry.mainfile}</TableCell>
              <TableCell>{new Date(entry.upload_time).toLocaleString()}</TableCell>
              <TableCell align="right">
                <DownloadButton query={{calc_id: entry.calc_id}} tooltip="Download files of this entry" />
                <Tooltip title="Show raw files and archive">
                  <IconButton onClick={() => history.push(`/entry/id/${entry.upload_id}/${entry.calc_id}`)}>
                    <MoreIcon />
                  </IconButton>
                </Tooltip>
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    )
  }
}

const Group = compose(withRouter, withApi(false), withStyles(GroupUnstyled.styles))(GroupUnstyled)

class GroupListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object,
    total: PropTypes.number,
    onChange: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired,
    groups_after: PropTypes.string,
    actions: PropTypes.element,
    domain: PropTypes.object.isRequired,
    selectedColumns: PropTypes.arrayOf(PropTypes.string)
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
    },
    details: {
      padding: 0
    }
  })

  addColumns(columns) {
    Object.keys(columns).forEach(key => {
      const column = columns[key]
      this.columns[key] = {
        ...column,
        supportsSort: false
      }
    })
  }

  constructor(props) {
    super(props)
    this.renderEntryActions = this.renderEntryActions.bind(this)

    this.columns = {}
  }

  componentDidMount() {
    this.addColumns(this.props.domain.searchResultColumns)
    this.addColumns(EntryListUnstyled.defaultColumns)
    this.addColumns({
      entries: {
        label: 'Entries',
        render: group => group.total.toLocaleString(),
        description: 'Number of entries in this group'
      }
    })
  }

  renderEntryActions(entry, selected) {
    return <DownloadButton
      dark={selected}
      query={{'dft.group_hash': entry.dft.group_hash}} tooltip="Download all entries of this group"
    />
  }

  renderEntryDetails(entry) {
    return <Group groupHash={entry.dft.group_hash} />
  }

  render() {
    const { classes, data, total, groups_after, onChange, actions, domain } = this.props
    const groups = data['dft.groups'] || {values: []}
    const results = Object.keys(groups.values).map(group_hash => {
      const example = groups.values[group_hash].examples[0]
      return {
        ...example,
        total: groups.values[group_hash].total,
        example: example
      }
    })
    const per_page = 10
    const after = groups.after

    const defaultSelectedColumns = this.props.selectedColumns || [
      ...domain.defaultSearchResultColumns,
      'datasets', 'authors', 'entries']

    let paginationText
    if (groups_after) {
      paginationText = `next ${results.length.toLocaleString()} of ${(total || 0).toLocaleString()}`
    } else {
      paginationText = `1-${results.length.toLocaleString()} of ${(total || 0).toLocaleString()}`
    }

    const pagination = <TableCell colSpan={1000} classes={{root: classes.scrollCell}}>
      <Toolbar className={classes.scrollBar}>
        <span className={classes.scrollSpacer}>&nbsp;</span>
        <span>{paginationText}</span>
        <IconButton disabled={!groups_after} onClick={() => onChange({groups_after: null})}>
          <StartIcon />
        </IconButton>
        <IconButton disabled={results.length < per_page} onClick={() => onChange({groups_after: after})}>
          <NextIcon />
        </IconButton>
      </Toolbar>
    </TableCell>

    return <DataTable
      classes={{details: classes.details}}
      entityLabels={['group of similar entries', 'groups of similar entries']}
      id={row => row.dft.group_hash}
      total={total}
      columns={this.columns}
      selectedColumns={defaultSelectedColumns}
      entryDetails={this.renderEntryDetails.bind(this)}
      entryActions={this.renderEntryActions}
      data={results}
      rows={per_page}
      actions={actions}
      pagination={pagination}
    />
  }
}

const GroupList = compose(withRouter, withApi(false), withStyles(GroupListUnstyled.styles))(GroupListUnstyled)

export default GroupList
