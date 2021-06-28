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
import FilterCheckboxes from './FilterCheckboxes'
import { useUnits } from '../../units'

const useFiltersElementStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  }
}))

export const labelElectronic = 'Electronic'
const options = {
  band_structure_electronic: {label: 'band structure'},
  dos_electronic: {label: 'density of states'}
}

const FilterElectronic = React.memo(({
  visible,
  className
}) => {
  const styles = useFiltersElementStyles()
  const units = useUnits()

  return <div className={clsx(className, styles.root)}>
    <Grid container spacing={2}>
      <Grid item xs={12}>
        <FilterSlider
          quantity="results.properties.electronic.band_structure_electronic.channel_info.band_gap"
          units={units}
          step={0.1}
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <FilterCheckboxes
          quantity="results.properties.electronic.band_structure_electronic.channel_info.band_gap_type"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <FilterCheckboxes
          label="available electronic properties"
          description="The electronic properties that are present in an entry."
          quantity="results.properties.available_properties"
          options={options}
          visible={visible}
        />
      </Grid>
    </Grid>
  </div>
})
FilterElectronic.propTypes = {
  visible: PropTypes.bool,
  className: PropTypes.string
}

export default FilterElectronic
