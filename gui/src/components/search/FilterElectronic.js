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
import React from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import { Grid } from '@material-ui/core'
import FilterSlider from './FilterSlider'
import { useUnits } from '../../units'

const useFiltersElementStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  }
}))

export const filterElectronic = [
  'results.properties.electronic.band_structure_electronic.channel_info.band_gap'
]
export const labelElectronic = 'Electronic'

/**
 * Displays the filter options for electronic properties.
 */
const FilterElectronic = React.memo(({
  className
}) => {
  const styles = useFiltersElementStyles()
  const units = useUnits()

  return <div className={clsx(className, styles.root)}>
    <Grid container spacing={2}>
      <Grid item xs={12}>
        <FilterSlider
          quantity="results.properties.electronic.band_structure_electronic.channel_info.band_gap"
          label="band gap"
          units={units}
        />
      </Grid>
    </Grid>
  </div>
})
FilterElectronic.propTypes = {
  className: PropTypes.string
}

export default FilterElectronic
