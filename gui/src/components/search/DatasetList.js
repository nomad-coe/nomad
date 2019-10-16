import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Table, TableHead, TableRow, TableCell, TableBody, Toolbar, IconButton } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import { withDomain } from '../domains'
import NextIcon from '@material-ui/icons/ChevronRight';
import StartIcon from '@material-ui/icons/SkipPrevious'


class DatasetListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    total: PropTypes.number.isRequired,
    onChange: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired,
    after: PropTypes.string
  }

  static styles = theme => ({
    root: {
      overflow: 'auto',
      paddingLeft: theme.spacing.unit * 2,
      paddingRight: theme.spacing.unit * 2,
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

  rowConfig = {
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
    },
  }

  handleClickDataset(dataset) {
    this.props.history.push(`/dataset/id/${dataset.id}`)
  }

  render() {
    const { classes, data, datasets_after, onChange } = this.props
    const results = Object.keys(data.datasets.values).map(id => ({
      id: id,
      total: data.datasets.values[id].total,
      example: data.datasets.values[id].examples[0]
    }))
    const per_page = 10
    const after = data.datasets.after

    const emptyRows = per_page - Math.min(per_page, results.length)

    return (
      <div className={classes.root}>
        <Table>
          <TableHead>
            <TableRow>
              {Object.keys(this.rowConfig).map(key => (
                <TableCell padding="dense" key={key}>
                  {this.rowConfig[key].label}
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {results.map((dataset, index) => (
              <TableRow hover tabIndex={-1} key={index} className={classes.clickableRow}>
                {Object.keys(this.rowConfig).map((key, rowIndex) => (
                  <TableCell padding="dense" key={rowIndex} onClick={() => this.handleClickDataset(dataset)} >
                    {this.rowConfig[key].render(dataset)}
                  </TableCell>
                ))}
              </TableRow>
            ))}
            {emptyRows > 0 && (
              <TableRow style={{ height: 57 * emptyRows }}>
                <TableCell colSpan={6} />
              </TableRow>
            )}
            <TableRow>
              <TableCell colSpan={1000} classes={{root: classes.scrollCell}}>
                <Toolbar className={classes.scrollBar}>
                  <span className={classes.scrollSpacer}>&nbsp;</span>
                  <IconButton disabled={!datasets_after} onClick={() => onChange({datasets_after: null})}>
                    <StartIcon />
                  </IconButton>
                  <IconButton onClick={() => onChange({datasets_after: after})}>
                    <NextIcon />
                  </IconButton>
                </Toolbar>
              </TableCell>
            </TableRow>
          </TableBody>
        </Table>
      </div>
    )
  }
}

const DatasetList = compose(withRouter, withDomain, withStyles(DatasetListUnstyled.styles))(DatasetListUnstyled)

Object.assign(DatasetList, {
  defaultState: {
    datasets_after: null
  }
})


export default DatasetList
