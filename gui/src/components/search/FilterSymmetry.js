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
import { Grid } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import FilterText from './FilterText'

export const labelSymmetry = 'Symmetry / Prototypes'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  }
}))

/**
 * Displays the filter options for symmetry properties.
 */
const FilterSymmetry = React.memo(({
  visible,
  className
}) => {
  const styles = useStyles()

  return <div className={clsx(className, styles.root)}>
    <Grid container spacing={2}>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.bravais_lattice"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.crystal_system"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.hall_symbol"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.point_group"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.space_group_symbol"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.prototype_aflow_id"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.structure_name"
          suggest
        />
      </Grid>
      <Grid item xs={12}>
        <FilterText
          quantity="results.material.symmetry.strukturbericht_designation"
          suggest
        />
      </Grid>
    </Grid>
  </div>
})
FilterSymmetry.propTypes = {
  visible: PropTypes.bool,
  className: PropTypes.string
}

export default FilterSymmetry
