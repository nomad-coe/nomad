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
  FormLabel, IconButton, MuiThemeProvider, Typography, Tooltip, TableSortLabel, ExpansionPanelDetails, ExpansionPanelSummary, ExpansionPanel, Grid } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from './errors'
import AnalyticsIcon from '@material-ui/icons/Settings'
import { analyticsTheme } from '../config'
import Link from 'react-router-dom/Link'
import { withApi } from './api'
import PeriodicTable from './PeriodicTable'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import QuantityHistogram from './QuantityHistogram'
import SearchBar from '../SearchBar'
import { withRouter } from 'react-router'

class Repo extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    raiseError: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    },
    searchDetails: {
      padding: 0,
      paddingBottom: theme.spacing.unit * 2,
      display: 'block',
      overflowX: 'auto'
    },
    searchSummary: {
      overflowX: 'auto'
    },
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
    },
    statistics: {
      minWidth: 500,
      maxWidth: 900,
      margin: 'auto',
      width: '100%'
    },
    quantityGrid: {
      minWidth: 524,
      maxWidth: 924,
      margin: 'auto',
      width: '100%'
    },
    quantity: {
      marginTop: theme.spacing.unit * 2
    },
    defaultCell: {
      overflow: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
      width: '100%',
      maxWidth: 200
    },
    searchBarContainer: {
      width: '100%',
      minWidth: 500,
      maxWidth: 900,
      margin: 'auto',
      marginBottom: theme.spacing.unit * 3
    }
  })

  rowConfig = {
    formula: {
      label: 'Formula'
    },
    code_name: {
      label: 'Code'
    },
    basis_set: {
      label: 'Basis set'
    },
    xc_functional: {
      label: 'XT treatment'
    },
    system: {
      label: 'System'
    },
    crystal_system: {
      label: 'Crystal system'
    },
    spacegroup_symbol: {
      label: 'Spacegroup'
    },
    authors: {
      label: 'Authors',
      render: (authors) => {
        if (authors.length > 3) {
          return authors.filter((_, index) => index < 2).map(author => author.name).join('; ') + 'et al'
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

  state = {
    data: [],
    page: 1,
    rowsPerPage: 10,
    total: 0,
    loading: true,
    owner: 'all',
    sortedBy: 'formula',
    sortOrder: 'asc',
    searchValues: {},
    aggregations: {},
    metrics: {},
    metric: 'code_runs'
  }

  update(changes) {
    changes = changes || {}
    const { page, rowsPerPage, owner, sortedBy, sortOrder, searchValues, metric } = {...this.state, ...changes}
    delete changes.metric
    this.setState({loading: true, ...changes})

    // code_runs is returned anyways
    const metrics_to_retrieve = metric === 'code_runs' ? [] : [metric]

    this.props.api.search({
      page: page,
      per_page: rowsPerPage,
      owner: owner || 'all',
      order_by: sortedBy,
      order: (sortOrder === 'asc') ? 1 : -1,
      total_metrics: metrics_to_retrieve,
      aggregation_metrics: metrics_to_retrieve,
      ...searchValues
    }).then(data => {
      const { pagination: { total, page, per_page }, results, aggregations, metrics } = data
      this.setState({
        data: results,
        aggregations: aggregations,
        metrics: metrics,
        page: page,
        rowsPerPage:
        per_page,
        total: total,
        loading: false,
        owner: owner,
        metric: metric
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
    this.update({page: page + 1})
  }

  handleChangeRowsPerPage = event => {
    const rowsPerPage = event.target.value
    this.update({rowsPerPage: rowsPerPage})
  }

  handleOwnerChange(owner) {
    this.update({owner: owner})
  }

  handleMetricChange(metric) {
    this.update({metric: metric})
  }

  handleSort(columnKey) {
    if (this.state.sortedBy === columnKey) {
      this.update({sortOrder: (this.state.sortOrder === 'asc') ? 'desc' : 'asc'})
    } else {
      this.update({sortOrder: 'asc', sortedBy: columnKey})
    }
  }

  handleClickCalc(calc) {
    this.props.history.push(`/repo/${calc.upload_id}/${calc.calc_id}`)
  }

  handleAtomsChanged(atoms) {
    const searchValues = {...this.state.searchValues}
    searchValues.atoms = atoms
    if (searchValues.atoms.length === 0) {
      delete searchValues.atoms
    }
    this.update({searchValues: searchValues})
  }

  handleQuantityChanged(quantity, selection) {
    const searchValues = {...this.state.searchValues}
    if (selection) {
      searchValues[quantity] = selection
    } else {
      delete searchValues[quantity]
    }
    this.update({searchValues: searchValues})
  }

  handleSearchChanged(searchValues) {
    this.update({searchValues: searchValues})
  }

  renderCell(key, rowConfig, calc) {
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

  render() {
    const { classes, user } = this.props
    const { data, rowsPerPage, page, total, loading, sortedBy, sortOrder, searchValues, aggregations, metrics, metric } = this.state
    const emptyRows = rowsPerPage - Math.min(rowsPerPage, total - (page - 1) * rowsPerPage)

    const quantity = (key, title) => (<QuantityHistogram
      classes={{root: classes.quantity}} title={title || key} width={300}
      data={aggregations[key]} metric={metric}
      value={searchValues[key]}
      onChanged={(selection) => this.handleQuantityChanged(key, selection)}/>)

    const ownerLabel = {
      all: 'All calculations',
      user: 'Your calculations',
      staging: 'Only calculations from your staging area'
    }

    const metricsLabel = {
      code_runs: 'Code runs',
      total_energies: 'Total energy calculations',
      geometries: 'Unique geometries',
      datasets: 'Datasets',
      unique_code_runs: 'Unique code runs'
    }
    return (
      <div className={classes.root}>
        { user
          ? <FormControl>
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
          </FormControl> : ''
        }

        <div className={classes.searchBarContainer}>
          <SearchBar
            fullWidth fullWidthInput={false} label="search" placeholder="enter atoms or other quantities"
            aggregations={aggregations} values={searchValues} metric={metric}
            onChanged={values => this.handleSearchChanged(values)}
          />
        </div>

        <ExpansionPanel>
          <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>} className={classes.searchSummary}>
            <Typography variant="h6" style={{textAlign: 'center', width: '100%', fontWeight: 'normal'}}>
              Found <b>{metrics.code_runs}</b>{metric === 'unique_code_runs' ? (<span>(<b>{metrics.unique_code_runs}</b> unique)</span>) : ''} code runs
              {metric === 'geometries' ? (<span> that simulate <b>{metrics.geometries}</b> unique geometries</span>) : ''}
              {metric === 'total_energies' ? (<span> with <b>{metrics.total_energies}</b> total energy calculations</span>) : ''}
              {metric === 'datasets' ? (<span> curated in <b>{metrics.datasets}</b> datasets</span>) : ''}.
            </Typography>
          </ExpansionPanelSummary>
          <ExpansionPanelDetails className={classes.searchDetails}>
            <div className={classes.statistics}>
              <FormControl>
                <FormLabel>Metric used in statistics: </FormLabel>
                <FormGroup row>
                  {['code_runs', 'unique_code_runs', 'total_energies', 'geometries', 'datasets'].map(metric => (
                    <FormControlLabel key={metric}
                      control={
                        <Checkbox checked={this.state.metric === metric} onChange={() => this.handleMetricChange(metric)} value={metric} />
                      }
                      label={metricsLabel[metric]}
                    />
                  ))}
                </FormGroup>
              </FormControl>
            </div>

            <PeriodicTable
              aggregations={aggregations.atoms} metric={metric}
              values={searchValues.atoms || []}
              onChanged={(selection) => this.handleAtomsChanged(selection)}
            />

            <Grid container spacing={24} className={classes.quantityGrid}>

              <Grid item xs={4}>
                {quantity('system', 'System')}
                {quantity('crystal_system', 'Crystal system')}
              </Grid>
              <Grid item xs={4}>
                {quantity('basis_set', 'Basis set')}
                {quantity('xc_functional', 'XC functionals')}
              </Grid>
              <Grid item xs={4}>
                {quantity('code_name', 'Code')}
              </Grid>
            </Grid>
          </ExpansionPanelDetails>
        </ExpansionPanel>

        <FormGroup className={classes.selectFormGroup} row>
          <FormLabel classes={{root: classes.selectLabel}} style={{flexGrow: 1}}>

          </FormLabel>
          <FormLabel classes={{root: classes.selectLabel}}>
            Analyse {total} code runs in an analytics notebook
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
                {Object.keys(this.rowConfig).map(key => (
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
                        {this.rowConfig[key].label}
                      </TableSortLabel>
                    </Tooltip>
                  </TableCell>
                ))}
              </TableRow>
            </TableHead>
            <TableBody>
              {data.map((calc, index) => (
                <TableRow hover tabIndex={-1} key={index} className={classes.clickableRow}>
                  {Object.keys(this.rowConfig).map((key, rowIndex) => (
                    <TableCell padding="dense" key={rowIndex} onClick={() => this.handleClickCalc(calc)} >
                      {this.renderCell(key, this.rowConfig[key], calc)}
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

export default compose(withRouter, withApi(false), withErrors, withStyles(Repo.styles))(Repo)
