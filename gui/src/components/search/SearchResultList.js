import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Paper, Table, TableHead, TableRow, TableCell, Tooltip, TableSortLabel, TableBody, TablePagination } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import { withDomain } from '../domains'

class SearchResultListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    total: PropTypes.number.isRequired,
    onChange: PropTypes.func,
    history: PropTypes.any.isRequired,
    order_by: PropTypes.string.isRequired,
    order: PropTypes.number.isRequired,
    page: PropTypes.number.isRequired,
    per_page: PropTypes.number.isRequired,
    domain: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      width: '100%',
      overflowX: 'scroll'
    },
    clickableRow: {
      cursor: 'pointer'
    }
  })

  defaultRowConfig = {
    authors: {
      label: 'Authors',
      render: (authors) => {
        if (authors.length > 3) {
          return authors.filter((_, index) => index < 2).map(author => author.name).join('; ') + ' et al'
        } else {
          return authors.map(author => author.name).join('; ')
        }
      }
    },
    references: {
      label: 'References',
      render: (references) => {
        if (references) {
          return references.map((reference, index) => (
            <div key={index} className={this.props.classes.defaultCell}>
              <a href={reference}>{reference}</a>
            </div>
          ))
        } else {
          return <i>no references</i>
        }
      }
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

  render() {
    const { classes, data, order, order_by, page, per_page, domain } = this.props
    const { results, pagination: { total } } = data

    const rowConfig = {
      ...domain.searchResultColumns,
      ...this.defaultRowConfig
    }

    const renderCell = (key, rowConfig, calc) => {
      const value = calc[key]
      if (rowConfig.render) {
        return rowConfig.render(value)
      } else {
        return (
          <div className={this.props.classes.defaultCell}>
            {value}
          </div>
        )
      }
    }

    const emptyRows = per_page - Math.min(per_page, total - (page - 1) * per_page)

    return (
      <Paper className={classes.root}>
        <Table>
          <TableHead>
            <TableRow>
              {Object.keys(rowConfig).map(key => (
                <TableCell padding="dense" key={key}>
                  <Tooltip
                    title="Sort"
                    placement={'bottom-start'}
                    enterDelay={300}
                  >
                    <TableSortLabel
                      active={order_by === key}
                      direction={order === 1 ? 'desc' : 'asc'}
                      onClick={() => this.handleSort(key)}
                    >
                      {rowConfig[key].label}
                    </TableSortLabel>
                  </Tooltip>
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {results.map((calc, index) => (
              <TableRow hover tabIndex={-1} key={index} className={classes.clickableRow}>
                {Object.keys(rowConfig).map((key, rowIndex) => (
                  <TableCell padding="dense" key={rowIndex} onClick={() => this.handleClickCalc(calc)} >
                    {renderCell(key, rowConfig[key], calc)}
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
              <TablePagination
                count={total}
                rowsPerPage={per_page}
                page={page - 1}
                backIconButtonProps={{
                  'aria-label': 'Previous Page'
                }}
                nextIconButtonProps={{
                  'aria-label': 'Next Page'
                }}
                onChangePage={this.handleChangePage}
                onChangeRowsPerPage={this.handleChangeRowsPerPage}
              />
            </TableRow>
          </TableBody>
        </Table>
      </Paper>
    )
  }
}

const SearchResultList = compose(withRouter, withDomain, withStyles(SearchResultListUnstyled.styles))(SearchResultListUnstyled)
Object.assign(SearchResultList, {
  defaultState: {
    order_by: 'formula',
    order: 1,
    page: 1,
    per_page: 10
  }
})

export default SearchResultList
