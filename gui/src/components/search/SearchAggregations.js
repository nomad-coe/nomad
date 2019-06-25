import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormControl, FormLabel, FormGroup, FormControlLabel, Checkbox, Tooltip } from '@material-ui/core'
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
    domain: PropTypes.object.isRequired,
    showDetails: PropTypes.bool
  }

  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit
    },
    details: {
      marginTop: theme.spacing.unit * 3
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
    const { classes, data, total_metrics, searchValues, domain, onChange, showDetails } = this.props
    const { aggregations, metrics } = data
    const selectedMetric = total_metrics.length === 0 ? 'code_runs' : total_metrics[0]
    const useMetric = Object.keys(metrics).find(metric => metric !== 'code_runs') || 'code_runs'
    const metricsDefinitions = domain.searchMetrics

    return (
      <div className={classes.root}>
        <div className={classes.details} style={showDetails ? {} : {display: 'none'}}>
          <FormControl>
            <FormLabel>Metric used in statistics: </FormLabel>
            <FormGroup row>
              {Object.keys(metricsDefinitions).map(metric => (
                <Tooltip key={metric} title={metricsDefinitions[metric].tooltip}>
                  <FormControlLabel
                    control={
                      <Checkbox checked={selectedMetric === metric} onChange={() => this.handleMetricChange(metric)} value={metric} />
                    }
                    label={metricsDefinitions[metric].label}
                  />
                </Tooltip>
              ))}
            </FormGroup>
          </FormControl>
          <domain.SearchAggregations aggregations={aggregations} searchValues={searchValues} metric={useMetric} onChange={onChange} />
        </div>
      </div>
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
