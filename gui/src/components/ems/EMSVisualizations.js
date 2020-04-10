import React from 'react'
import { Grid } from '@material-ui/core'
import { Quantity } from '../search/QuantityHistogram'
import SearchContext from '../search/SearchContext'

export default class EMSVisualizations extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {usedMetric}} = this.context

    return (
      <Grid container spacing={24}>
        <Grid item xs={6}>
          <Quantity quantity="ems.method" title="Method" scale={1} metric={usedMetric} />
          <Quantity quantity="ems.probing_method" title="Probing" scale={1} metric={usedMetric} />
        </Grid>
        <Grid item xs={6}>
          <Quantity quantity="ems.sample_microstructure" title="Sample structure" scale={1} metric={usedMetric} />
          <Quantity quantity="ems.sample_constituents" title="Sample constituents" scale={1} metric={usedMetric} />
        </Grid>
      </Grid>
    )
  }
}
