import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormControl, FormLabel, FormGroup, FormControlLabel, Checkbox, Tooltip } from '@material-ui/core'
import { withDomain } from '../domains'
import { compose } from 'recompose'
import SearchContext from '../search/SearchContext'

class SearchAggregations extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
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

  static contextType = SearchContext.type

  render() {
    const {classes, domain, showDetails} = this.props
    const {state: {response: {statistics}, metric}, setMetric} = this.context
    const selectedMetric = metric

    // first the first statistic to determine which metric is used
    let useMetric = 'code_runs'
    const firstRealQuantitiy = Object.keys(statistics).find(key => key !== 'total')
    if (firstRealQuantitiy) {
      const firstValue = Object.keys(statistics[firstRealQuantitiy])[0]
      if (firstValue) {
        useMetric = Object.keys(statistics[firstRealQuantitiy][firstValue])
          .find(metric => metric !== 'code_runs') || 'code_runs'
      }
    }

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
                      <Checkbox checked={selectedMetric === metric} onChange={() => setMetric(metric)} value={metric} />
                    }
                    label={metricsDefinitions[metric].label}
                  />
                </Tooltip>
              ))}
            </FormGroup>
          </FormControl>
          <domain.SearchAggregations metric={useMetric} />
        </div>
      </div>
    )
  }
}

export default compose(withDomain, withStyles(SearchAggregations.styles))(SearchAggregations)
