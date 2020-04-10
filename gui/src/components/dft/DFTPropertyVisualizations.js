import React from 'react'
import PropTypes from 'prop-types'
import { Grid } from '@material-ui/core'
import { Quantity } from '../search/QuantityHistogram'
import SearchContext from '../search/SearchContext'

export default class DFTPropertyVisualizations extends React.Component {
  static propTypes = {
    info: PropTypes.object
  }

  static contextType = SearchContext.type

  componentDidMount() {
    const {setStatisticsToRefresh} = this.context
    setStatisticsToRefresh('dft.labels_springer_classification')
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
          <Quantity quantity="dft.quantities_energy" title="Energy" scale={1} metric={usedMetric} sort tooltips />
          <Quantity quantity="dft.quantities_electronic" title="Electronic" scale={1} metric={usedMetric} sort tooltips />
        </Grid>
        <Grid item xs={4}>
          <Quantity quantity="dft.quantities_forces" title="Forces" scale={1} metric={usedMetric} sort tooltips />
          <Quantity quantity="dft.quantities_vibrational" title="Vibrational" scale={1} metric={usedMetric} sort tooltips />
          <Quantity quantity="dft.quantities_optical" title="Optical" scale={1} metric={usedMetric} sort tooltips />
        </Grid>
        <Grid item xs={4}>
          <Quantity quantity="dft.labels_springer_classification" title="Springer classification" scale={1} metric={usedMetric} tooltips />
          <Quantity quantity="dft.quantities_magnetic" title="Magnetic" scale={1} metric={usedMetric} sort tooltips />
        </Grid>
      </Grid>
    )
  }
}
