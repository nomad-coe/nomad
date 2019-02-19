import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Table from '@material-ui/core/Table'
import TableBody from '@material-ui/core/TableBody'
import TableCell from '@material-ui/core/TableCell'
import TablePagination from '@material-ui/core/TablePagination'
import TableRow from '@material-ui/core/TableRow'
import Paper from '@material-ui/core/Paper'
import { TableHead, LinearProgress, FormControl, FormControlLabel, Checkbox, FormGroup,
  FormLabel, IconButton, MuiThemeProvider, Typography, Tooltip, TableSortLabel } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from './errors'
import AnalyticsIcon from '@material-ui/icons/Settings'
import { analyticsTheme } from '../config'
import Link from 'react-router-dom/Link'
import { withApi } from './api'
import CalcDialog from './CalcDialog'
// import PeriodicTable from './PeriodicTable'

class Repo extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    data: {
      width: '100%',
      overflowX: 'scroll'
    },
    title: {
      marginBottom: theme.spacing.unit * 4
    },
    progressPlaceholder: {
      height: 5
    },
    summary: {
      textAlign: 'center',
      marginTop: theme.spacing.unit * 2
    },
    selectFormGroup: {
      paddingLeft: theme.spacing.unit * 3
    },
    selectLabel: {
      padding: theme.spacing.unit * 2
    },
    clickableRow: {
      cursor: 'pointer'
    }
  })

  static rowConfig = {
    formula: 'Formula',
    code_name: 'Code',
    basis_set: 'Basis set',
    system: 'System',
    crystal_system: 'Crystal system',
    spacegroup: 'Spacegroup',
    xc_functional: 'XT treatment'
  }

  state = {
    data: [],
    page: 1,
    rowsPerPage: 10,
    total: 0,
    loading: true,
    owner: 'all',
    sortedBy: 'formula',
    sortOrder: 'asc',
    openCalc: null
  }

  update(changes) {
    changes = changes || {}
    const { page, rowsPerPage, owner, sortedBy, sortOrder } = {...this.state, ...changes}
    this.setState({loading: true, ...changes})

    this.props.api.search({
      page: page,
      per_page: rowsPerPage,
      owner: owner || 'all',
      order_by: sortedBy,
      order: (sortOrder === 'asc') ? 1 : -1
    }).then(data => {
      const { pagination: { total, page, per_page }, results } = data
      this.setState({
        data: results,
        page: page,
        rowsPerPage:
        per_page,
        total: total,
        loading: false,
        owner: owner
      })
    }).catch(errors => {
      this.setState({data: [], total: 0, loading: false, owner: owner})
      this.props.raiseError(errors)
    })
  }

  componentDidMount() {
    this.update()
  }

  handleChangePage = (event, page) => {
    this.update({page: this.state.page + 1})
  }

  handleChangeRowsPerPage = event => {
    const rowsPerPage = event.target.value
    this.update({rowsPerPage: rowsPerPage})
  }

  handleOwnerChange(owner) {
    this.update({owner: owner})
  }

  handleSort(columnKey) {
    if (this.state.sortedBy === columnKey) {
      this.update({sortOrder: (this.state.sortOrder === 'asc') ? 'desc' : 'asc'})
    } else {
      this.update({sortOrder: 'asc', sortedBy: columnKey})
    }
  }

  handleCalcClose() {
    this.setState({openCalc: null})
  }

  handleClickCalc(calc_id) {
    this.setState({openCalc: calc_id})
  }

  render() {
    const { classes } = this.props
    const { data, rowsPerPage, page, total, loading, sortedBy, sortOrder, openCalc } = this.state
    const emptyRows = rowsPerPage - Math.min(rowsPerPage, total - (page - 1) * rowsPerPage)

    const ownerLabel = {
      all: 'All calculations',
      user: 'Your calculations',
      staging: 'Only calculations from your staging area'
    }
    return (
      <div className={classes.root}>
        { openCalc ? <CalcDialog calcId={openCalc.calc_id} uploadId={openCalc.upload_id} onClose={() => this.handleCalcClose()} /> : ''}
        <Typography variant="h4" className={classes.title}>The Repository – Raw Code Data</Typography>
        {/* <PeriodicTable/> */}
        <FormControl>
          <FormLabel>Filter calculations and only show: </FormLabel>
          <FormGroup row>
            {['all', 'user', 'staging'].map(owner => (
              <FormControlLabel key={owner}
                control={
                  <Checkbox checked={this.state.owner === owner} onChange={() => this.handleOwnerChange(owner)} value="owner" />
                }
                label={ownerLabel[owner]}
              />
            ))}
          </FormGroup>
        </FormControl>

        <FormGroup className={classes.selectFormGroup} row>
          <FormLabel classes={{root: classes.selectLabel}} style={{flexGrow: 1}}>

          </FormLabel>
          <FormLabel classes={{root: classes.selectLabel}}>
            Analyse {total} calculations in an analytics notebook
          </FormLabel>
          <MuiThemeProvider theme={analyticsTheme}>
            <IconButton color="primary" component={Link} to={`/analytics`}>
              <AnalyticsIcon />
            </IconButton>
          </MuiThemeProvider>
        </FormGroup>

        <Paper className={classes.data}>
          {loading ? <LinearProgress variant="query" /> : <div className={classes.progressPlaceholder} />}
          <Table>
            <TableHead>
              <TableRow>
                {Object.keys(Repo.rowConfig).map(key => (
                  <TableCell padding="dense" key={key}>
                    <Tooltip
                      title="Sort"
                      placement={'bottom-start'}
                      enterDelay={300}
                    >
                      <TableSortLabel
                        active={sortedBy === key}
                        direction={sortOrder}
                        onClick={() => this.handleSort(key)}
                      >
                        {Repo.rowConfig[key]}
                      </TableSortLabel>
                    </Tooltip>
                  </TableCell>
                ))}
              </TableRow>
            </TableHead>
            <TableBody>
              {data.map((calc, index) => (
                <TableRow hover tabIndex={-1} key={index} className={classes.clickableRow}>
                  {Object.keys(Repo.rowConfig).map((key, rowIndex) => (
                    <TableCell padding="dense" key={rowIndex} onClick={() => this.handleClickCalc(calc)} >
                      {calc[key]}
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
                  rowsPerPage={rowsPerPage}
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
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(Repo.styles))(Repo)
