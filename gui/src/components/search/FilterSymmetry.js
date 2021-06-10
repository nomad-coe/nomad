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
import { useUnits } from '../../units'
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
  className
}) => {
  const styles = useStyles()
  const units = useUnits()

  return <div className={clsx(className, styles.root)}>
    <Grid container spacing={2}>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.bravais_lattice"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.crystal_system"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.hall_symbol"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.point_group"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.space_group_symbol"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.prototype_aflow_id"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.structure_name"
          units={units}
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.symmetry.strukturbericht_designation"
          units={units}
        />
      </Grid>
    </Grid>
  </div>
})
FilterSymmetry.propTypes = {
  className: PropTypes.string
}

export default FilterSymmetry
