import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, ExpansionPanel, ExpansionPanelSummary, ExpansionPanelDetails,
  FormControl, FormLabel, FormGroup, FormControlLabel, Checkbox, Typography
} from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import { withDomain } from '../domains'
import { compose } from 'recompose'

class SearchAggregationsUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    onChange: PropTypes.func.isRequired,
    data: PropTypes.object.isRequired,
    total_metrics: PropTypes.arrayOf(PropTypes.string).isRequired,
    aggregation_metrics: PropTypes.arrayOf(PropTypes.string).isRequired,
    searchValues: PropTypes.object.isRequired,
    domain: PropTypes.object.isRequired
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
    }
  })

  handleMetricChange(metric) {
    const metrics = metric === 'code_runs' ? [] : [metric]
    this.setState({metric: metric})
    this.props.onChange({total_metrics: metrics, aggregation_metrics: metrics})
  }

  handleSearchChanged(searchValues) {
    this.props.onChange({searchValues: searchValues})
  }

  render() {
    const { classes, data, total_metrics, searchValues, domain, onChange } = this.props
    const { aggregations, metrics } = data
    const selectedMetric = total_metrics.length === 0 ? 'code_runs' : total_metrics[0]
    const useMetric = Object.keys(metrics).find(metric => metric !== 'code_runs') || 'code_runs'

    const metricsDefinitions = {
      code_runs: {
        label: 'Entries',
        renderResultString: count => (<span><b>{count}</b> entries</span>)
      },
      unique_code_runs: {
        label: 'Unique entries',
        renderResultString: count => (<span> and <b>{count}</b> unique entries</span>)
      },
      ...domain.searchMetrics,
      datasets: {
        label: 'Datasets',
        renderResultString: count => (<span> curated in <b>{metrics.datasets}</b> datasets</span>)
      }
    }

    return (
      <ExpansionPanel className={classes.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>} className={classes.searchSummary}>
          <Typography variant="h6" style={{textAlign: 'center', width: '100%', fontWeight: 'normal'}}>
            Found {Object.keys(metricsDefinitions).map(key => {
              return (key === useMetric || key === 'code_runs') ? <span key={key}>{metricsDefinitions[key].renderResultString(metrics[key])}</span> : ''
            })}.
          </Typography>
        </ExpansionPanelSummary>
        <ExpansionPanelDetails className={classes.searchDetails}>
          <div className={classes.statistics}>
            <FormControl>
              <FormLabel>Metric used in statistics: </FormLabel>
              <FormGroup row>
                {Object.keys(metricsDefinitions).map(metric => (
                  <FormControlLabel key={metric}
                    control={
                      <Checkbox checked={selectedMetric === metric} onChange={() => this.handleMetricChange(metric)} value={metric} />
                    }
                    label={metricsDefinitions[metric].label}
                  />
                ))}
              </FormGroup>
            </FormControl>
          </div>

          <domain.SearchAggregations aggregations={aggregations} searchValues={searchValues} metric={useMetric} onChange={onChange} />
        </ExpansionPanelDetails>
      </ExpansionPanel>
    )
  }
}

const SearchAggregations = compose(withDomain, withStyles(SearchAggregationsUnstyled.styles))(SearchAggregationsUnstyled)
Object.assign(SearchAggregations, {
  defaultState: {
    aggregation_metrics: [],
    total_metrics: [],
    searchValues: {}
  }
})

export default SearchAggregations
