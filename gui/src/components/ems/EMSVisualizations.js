/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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
