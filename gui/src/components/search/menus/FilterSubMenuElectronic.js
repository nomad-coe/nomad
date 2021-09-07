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
import React, { useContext } from 'react'
import PropTypes from 'prop-types'
import { Grid } from '@material-ui/core'
import { FilterSubMenu, filterMenuContext } from './FilterMenu'
import InputSlider from '../input/InputSlider'
import InputCheckboxes from '../input/InputCheckboxes'
import { Quantity, useUnits } from '../../../units'

const step = new Quantity(0.1, 'electron_volt')
const options = {
  band_structure_electronic: {label: 'band structure'},
  dos_electronic: {label: 'density of states'}
}

const FilterSubMenuElectronic = React.memo(({
  value,
  ...rest
}) => {
  const units = useUnits()
  const {selected} = useContext(filterMenuContext)
  const visible = value === selected

  return <FilterSubMenu value={value} {...rest}>
    <Grid container spacing={2}>
      <Grid item xs={12}>
        <InputSlider
          quantity="results.properties.electronic.band_structure_electronic.channel_info.band_gap"
          units={units}
          step={step}
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <InputCheckboxes
          quantity="results.properties.electronic.band_structure_electronic.channel_info.band_gap_type"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <InputCheckboxes
          label="available electronic properties"
          description="The electronic properties that are present in an entry."
          quantity="electronic_properties"
          options={options}
          visible={visible}
        />
      </Grid>
    </Grid>
  </FilterSubMenu>
})
FilterSubMenuElectronic.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuElectronic
