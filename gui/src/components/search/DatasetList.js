import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, TableCell, Toolbar, IconButton } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import { withDomain } from '../domains'
import NextIcon from '@material-ui/icons/ChevronRight'
import StartIcon from '@material-ui/icons/SkipPrevious'
import DataTable from '../DataTable'

class DatasetListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    total: PropTypes.number.isRequired,
    onChange: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired,
    datasets_after: PropTypes.string
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

  columns = {
    name: {
      label: 'Dataset name',
      render: (dataset) => dataset.example.datasets.find(d => d.id + '' === dataset.id).name
    },
    DOI: {
      label: 'Dataset DOI',
      render: (dataset) => dataset.example.datasets.find(d => d.id + '' === dataset.id).doi
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

  handleClickDataset(datasetId) {
    this.props.history.push(`/dataset/id/${datasetId}`)
  }

  render() {
    const { classes, data, total, datasets_after, onChange } = this.props
    const results = Object.keys(data.datasets.values).map(id => ({
      id: id,
      total: data.datasets.values[id].total,
      example: data.datasets.values[id].examples[0]
    }))
    const per_page = 10
    const after = data.datasets.after

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
      // entryActions={this.renderEntryActions.bind(this)}
      onEntryClicked={this.handleClickDataset.bind(this)}
      data={results}
      rows={per_page}
      pagination={pagination}
    />
  }
}

const DatasetList = compose(withRouter, withDomain, withStyles(DatasetListUnstyled.styles))(DatasetListUnstyled)

Object.assign(DatasetList, {
  defaultState: {
    datasets_after: null
  }
})

export default DatasetList
