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
import FilterCheckboxes from './FilterCheckboxes'

const useFiltersElementStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  }
}))

export const labelDFT = 'DFT'

/**
 * Displays the filter options for electronic properties.
 */
const FilterDFT = React.memo(({
  visible,
  className
}) => {
  const styles = useFiltersElementStyles()

  return <div className={clsx(className, styles.root)}>
    <Grid container spacing={2}>
      <Grid item xs={12}>
        <FilterCheckboxes
          quantity="results.method.simulation.dft.basis_set_type"
          visible={visible}
          xs={6}
        />
      </Grid>
      <Grid item xs={12}>
        <FilterCheckboxes
          quantity="results.method.simulation.dft.core_electron_treatment"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <FilterCheckboxes
          quantity="results.method.simulation.dft.relativity_method"
          visible={visible}
        />
      </Grid>
      {/* <Grid item xs={12}>
        <FilterSelect
          quantity="results.method.simulation.dft.basis_set_name"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <FilterCheckboxes
          quantity="results.method.simulation.dft.van_der_Waals_method"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <FilterCheckboxes
          quantity="results.method.simulation.dft.smearing_type"
          visible={visible}
        />
      </Grid> */}
    </Grid>
  </div>
})
FilterDFT.propTypes = {
  visible: PropTypes.bool,
  className: PropTypes.string
}

export default FilterDFT
