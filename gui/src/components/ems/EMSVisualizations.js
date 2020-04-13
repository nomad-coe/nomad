import React, { useContext, useEffect } from 'react'
import { Grid } from '@material-ui/core'
import { Quantity } from '../search/QuantityHistogram'
import { searchContext } from '../search/SearchContext'

export default function EMSVisualizations(props) {
  const {state: {usedMetric}, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['ems.method', 'ems.probing_method', 'ems.sample_microstructure', 'ems.sample_constituents'])
  }, [])
  return (
    <Grid container spacing={2}>
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
