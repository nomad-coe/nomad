import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Grid } from '@material-ui/core'
import QuantityHistogram from '../search/QuantityHistogram'
import SearchContext from '../search/SearchContext'
import { withApi } from '../api'


class QuantityUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    quantity: PropTypes.string.isRequired,
    metric: PropTypes.string.isRequired,
    title: PropTypes.string,
    scale: PropTypes.number
  }
  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2
    }
  })

  static contextType = SearchContext.type

  render() {
    const {classes, scale, quantity, title, ...props} = this.props
    const {state: {response, query}, setQuery} = this.context

    return <QuantityHistogram
      classes={{root: classes.root}}
      width={300}
      defaultScale={scale || 1}
      title={title || quantity}
      data={response.statistics[quantity]}
      value={query[quantity]}
      onChanged={selection => setQuery({...query, [quantity]: selection})}
      {...props} />
  }
}

const Quantity = withStyles(QuantityUnstyled.styles)(QuantityUnstyled)


class DFTSearchAggregations extends React.Component {
  static propTypes = {
    info: PropTypes.object
  }

  static contextType = SearchContext.type

  render() {
    const {info} = this.props
    const {state: {response: {statistics}, usedMetric}} = this.context

    if (statistics.code_name && info) {
      // filter based on known codes, since elastic search might return 0 aggregations on
      // obsolete code names
      const filteredCodeNames = {}
      const defaultValue = {
        code_runs: 0
      }
      defaultValue[usedMetric] = 0
      info.codes.forEach(key => {
        filteredCodeNames[key] = statistics.code_name[key] || defaultValue
      })
      statistics.code_name = filteredCodeNames
    }

    return (
      <Grid container spacing={24}>
        <Grid item xs={4}>
          <Quantity quantity="code_name" title="Code" scale={0.25} metric={usedMetric} />
        </Grid>
        <Grid item xs={4}>
          <Quantity quantity="system" title="System type" scale={0.25} metric={usedMetric} />
          <Quantity quantity="crystal_system" title="Crystal system" scale={1} metric={usedMetric} />
        </Grid>
        <Grid item xs={4}>
          <Quantity quantity="basis_set" title="Basis set" scale={0.25} metric={usedMetric} />
          <Quantity quantity="xc_functional" title="XC functionals" scale={0.5} metric={usedMetric} />
        </Grid>
      </Grid>
    )
  }
}

export default withApi(false, false)(DFTSearchAggregations)
