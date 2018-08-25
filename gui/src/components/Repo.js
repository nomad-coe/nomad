import React from 'react';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TablePagination from '@material-ui/core/TablePagination';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import { TableFooter } from '@material-ui/core';
import api from '../api';


const styles = theme => ({
  root: {
    width: '100%',
    marginTop: theme.spacing.unit * 3,
  },
  table: {
  },
  tableWrapper: {
    // overflowX: 'auto',
  },
});

class EnhancedTable extends React.Component {
  state = {
    data: [],
    page: 1,
    rowsPerPage: 5,
    total: 0
  };

  update(page, rowsPerPage) {
    api.repoAll(page, rowsPerPage).then(data => {
      const { pagination: { total, page, per_page }, results } = data
      this.setState({data: results, page: page, rowsPerPage: per_page, total: total})
    })
  }

  componentDidMount() {
    const { page, rowsPerPage } = this.state
    this.update(page, rowsPerPage)
  }

  handleChangePage = (event, page) => {
    this.update(page + 1, this.state.rowsPerPage)
  };

  handleChangeRowsPerPage = event => {
    const rowsPerPage = event.target.value
    this.update(this.state.page, rowsPerPage)
  };

  render() {
    const { classes } = this.props;
    const { data, rowsPerPage, page, total } = this.state;
    const emptyRows = rowsPerPage - Math.min(rowsPerPage, total - (page - 1) * rowsPerPage);

    return (
      <Paper className={classes.root}>
        <div className={classes.tableWrapper}>
          <Table className={classes.table} aria-labelledby="tableTitle">
            <TableBody>
              {data.map((n, index) => (
                <TableRow
                  hover
                  role="checkbox"
                  tabIndex={-1}
                  key={index}
                >
                  <TableCell>{n.chemical_composition_bulk_reduced}</TableCell>
                  <TableCell>{n.program_name}</TableCell>
                  <TableCell>{n.program_basis_set_type}</TableCell>
                  <TableCell>{n.system_type}</TableCell>
                  <TableCell>{n.crystal_system}</TableCell>
                  <TableCell>{n.space_group_number}</TableCell>
                  <TableCell>{n.XC_functional_name}</TableCell>
                </TableRow>
              ))}
              {emptyRows > 0 && (
                <TableRow style={{ height: 49 * emptyRows }}>
                  <TableCell colSpan={6} />
                </TableRow>
              )}
            </TableBody>
            <TableFooter>
              <TablePagination
                count={total}
                rowsPerPage={rowsPerPage}
                page={page - 1}
                backIconButtonProps={{
                  'aria-label': 'Previous Page',
                }}
                nextIconButtonProps={{
                  'aria-label': 'Next Page',
                }}
                onChangePage={this.handleChangePage}
                onChangeRowsPerPage={this.handleChangeRowsPerPage}
              />
            </TableFooter>
          </Table>
        </div>
      </Paper>
    );
  }
}

EnhancedTable.propTypes = {
  classes: PropTypes.object.isRequired,
};

export default withStyles(styles)(EnhancedTable);