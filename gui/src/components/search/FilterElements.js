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
import React, { useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import { Grid } from '@material-ui/core'
import NewPeriodicTable from './NewPeriodicTable'
import FilterText from './FilterText'
import FilterSlider from './FilterSlider'
import { useFilterState, useAgg, useExclusiveState } from './FilterContext'
import { useUnits } from '../../units'

export const labelElements = 'Elements / Formula'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    marginTop: theme.spacing(0.5)
  },
  grid: {
    marginTop: theme.spacing(2)
  }
}))

const FilterElements = React.memo(({
  visible,
  className
}) => {
  const styles = useStyles()
  const [exclusive, setExclusive] = useExclusiveState()
  const [filter, setFilter] = useFilterState('results.material.elements')
  const data = useAgg('results.material.elements', 'terms', false, visible)
  const units = useUnits()
  const availableValues = useMemo(() => {
    const elementCountMap = {}
    if (data) {
      for (let value of data) {
        elementCountMap[value.value] = value.count
      }
    }
    return elementCountMap
  }, [data])

  const handleElementsChanged = useCallback(elements => {
    setFilter(elements)
  }, [setFilter])

  // If this panel is not visible, we hide the periodic table as it is very
  // expensive to re-render (it is rerendered always when the query changes).
  return <div className={clsx(className, styles.root)}>
    {visible &&
    <NewPeriodicTable
      availableValues={availableValues}
      values={filter}
      exclusive={exclusive}
      onChanged={handleElementsChanged}
      onExclusiveChanged={setExclusive}
    />}
    <Grid container spacing={2} className={styles.grid}>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.chemical_formula_hill"
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.chemical_formula_anonymous"
        />
      </Grid>
      <Grid item xs={12}>
        <FilterSlider
          quantity="results.material.n_elements"
          step={1}
          units={units}
          visible={visible}
        />
      </Grid>
    </Grid>
  </div>
})
FilterElements.propTypes = {
  visible: PropTypes.bool,
  className: PropTypes.string
}

export default FilterElements
