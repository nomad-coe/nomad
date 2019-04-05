import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, ExpansionPanel, ExpansionPanelSummary, ExpansionPanelDetails,
  FormControl, FormLabel, FormGroup, FormControlLabel, Checkbox, Grid, Typography
} from '@material-ui/core'
import PeriodicTable from './PeriodicTable'
import QuantityHistogram from './QuantityHistogram'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'

class SearchStatistics extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    onChange: PropTypes.func,
    data: PropTypes.object.isRequired,
    total_metrics: PropTypes.arrayOf(PropTypes.string).isRequired,
    aggregation_metrics: PropTypes.arrayOf(PropTypes.string).isRequired,
    searchValues: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {},
    searchDetails: {
      padding: 0,
      paddingBottom: theme.spacing.unit * 2,
      display: 'block',
      overflowX: 'auto'
    },
    searchSummary: {
      overflowX: 'auto'
    },
    summary: {
      textAlign: 'center',
      marginTop: theme.spacing.unit * 2
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
    }
  })

  static defaultState = ({
    aggregation_metrics: [],
    total_metrics: [],
    searchValues: {}
  })

  metricsLabels = {
    code_runs: 'Entries',
    total_energies: 'Total energy calculations',
    geometries: 'Unique geometries',
    datasets: 'Datasets',
    unique_code_runs: 'Unique entries'
  }

  handleChange(changes) {
    if (this.props.onChange) {
      this.props.onChange(changes)
    }
  }

  handleMetricChange(metric) {
    const metrics = metric === 'code_runs' ? [] : [metric]
    this.setState({metric: metric})
    this.handleChange({total_metrics: metrics, aggregation_metrics: metrics})
  }

  handleAtomsChanged(atoms) {
    const searchValues = {...this.props.searchValues}
    searchValues.atoms = atoms
    if (searchValues.atoms.length === 0) {
      delete searchValues.atoms
    }
    this.handleChange({searchValues: searchValues})
  }

  handleQuantityChanged(quantity, selection) {
    const searchValues = {...this.props.searchValues}
    if (selection) {
      searchValues[quantity] = selection
    } else {
      delete searchValues[quantity]
    }
    this.handleChange({searchValues: searchValues})
  }

  handleSearchChanged(searchValues) {
    this.handleChange({searchValues: searchValues})
  }

  render() {
    const { classes, data, total_metrics, searchValues } = this.props
    const { aggregations, metrics } = data
    const selectedMetric = total_metrics.length === 0 ? 'code_runs' : total_metrics[0]
    const useMetric = Object.keys(metrics).find(metric => metric !== 'code_runs') || 'code_runs'

    const quantity = (key, title) => (<QuantityHistogram
      classes={{root: classes.quantity}} title={title || key} width={300}
      data={aggregations[key]} metric={useMetric}
      value={searchValues[key]}
      onChanged={(selection) => this.handleQuantityChanged(key, selection)}/>)

    return (
      <ExpansionPanel className={classes.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>} className={classes.searchSummary}>
          <Typography variant="h6" style={{textAlign: 'center', width: '100%', fontWeight: 'normal'}}>
            Found <b>{metrics.code_runs}</b>{useMetric === 'unique_code_runs' ? (<span>(<b>{metrics.unique_code_runs}</b> unique)</span>) : ''} code runs
            {useMetric === 'geometries' ? (<span> that simulate <b>{metrics.geometries}</b> unique geometries</span>) : ''}
            {useMetric === 'total_energies' ? (<span> with <b>{metrics.total_energies}</b> total energy calculations</span>) : ''}
            {useMetric === 'datasets' ? (<span> curated in <b>{metrics.datasets}</b> datasets</span>) : ''}.
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
                      <Checkbox checked={selectedMetric === metric} onChange={() => this.handleMetricChange(metric)} value={metric} />
                    }
                    label={this.metricsLabels[metric]}
                  />
                ))}
              </FormGroup>
            </FormControl>
          </div>

          <PeriodicTable
            aggregations={aggregations.atoms} metric={useMetric}
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
    )
  }
}

export default withStyles(SearchStatistics.styles)(SearchStatistics)
