import React, { useContext, useEffect } from 'react'
import { Grid } from '@material-ui/core'
import QuantityHistogram from '../search/QuantityHistogram'
import { searchContext } from '../search/SearchContext'

export default function EMSVisualizations(props) {
  const {setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['ems.method', 'ems.probing_method', 'ems.sample_microstructure', 'ems.sample_constituents'])
    // eslint-disable-next-line
  }, [])

  return (
    <Grid container spacing={2}>
      <Grid item xs={6}>
        <QuantityHistogram quantity="ems.method" title="Method" />
        <QuantityHistogram quantity="ems.probing_method" title="Probing" />
      </Grid>
      <Grid item xs={6}>
        <QuantityHistogram quantity="ems.sample_microstructure" title="Sample structure" />
        <QuantityHistogram quantity="ems.sample_constituents" title="Sample constituents" />
      </Grid>
    </Grid>
  )
}
