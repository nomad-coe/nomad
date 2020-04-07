import React from 'react'
import PropTypes from 'prop-types'
import { Grid } from '@material-ui/core'
import { Quantity } from '../search/QuantityHistogram'
import SearchContext from '../search/SearchContext'
import { withApi } from '../api'

class DFTSearchAggregations extends React.Component {
  static propTypes = {
    info: PropTypes.object
  }

  static contextType = SearchContext.type

  componentDidMount() {
    const {setStatisticsToRefresh} = this.context
    setStatisticsToRefresh('dft.labels_springer_compound_class')
  }

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
          <Quantity quantity="dft.code_name" title="Code" scale={0.25} metric={usedMetric} />
        </Grid>
        <Grid item xs={4}>
          <Quantity quantity="dft.basis_set" title="Basis set" scale={0.25} metric={usedMetric} />
          <Quantity quantity="dft.xc_functional" title="XC functionals" scale={0.5} metric={usedMetric} />
          <Quantity quantity="dft.compound_type" title="Compound type" scale={1} metric={usedMetric} />
        </Grid>
        <Grid item xs={4}>
          <Quantity quantity="dft.system" title="System type" scale={0.25} metric={usedMetric} />
          <Quantity quantity="dft.crystal_system" title="Crystal system" scale={1} metric={usedMetric} />
          <Quantity quantity="dft.labels_springer_compound_class" title="Springer compound class" scale={1} metric={usedMetric} />
        </Grid>
      </Grid>
    )
  }
}

export default withApi(false, false)(DFTSearchAggregations)
