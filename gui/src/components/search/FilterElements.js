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
import React, { useState, useContext, useEffect, useCallback } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import { Grid } from '@material-ui/core'
import NewPeriodicTable from './NewPeriodicTable'
import FilterText from './FilterText'
import { searchContext } from './SearchContext'
import { useFilters } from './FilterContext'

const useFiltersElementStyles = makeStyles(theme => ({
  root: {
    width: '40rem',
    marginTop: theme.spacing(0.5)
  },
  grid: {
    marginTop: theme.spacing(2)
  }

}))

/**
 * Displays the filter options for chemical elements.
 */
const FilterElements = React.memo(({
  className
}) => {
  const styles = useFiltersElementStyles()
  const [exclusive, setExclusive] = useState(false)
  const {response: {statistics, metric}, query, setQuery, setStatistics} = useContext(searchContext)
  const {elements, setElements} = useFilters()

  useEffect(() => {
    setStatistics(['atoms'])
    // eslint-disable-next-line
  }, [])

  const handleExclusiveChanged = () => {
    // const newExclusive = !exclusive
    // if (newExclusive) {
    //   setQuery({only_atoms: query.atoms, atoms: []})
    // } else {
    //   setQuery({atoms: query.only_atoms, only_atoms: []})
    // }
    // setExclusive(newExclusive)
  }

  const handleElementsChanged = useCallback(elements => {
    setElements(new Set(elements))
    // exclusive && setExclusive(false)
    // setQuery({atoms: atoms, only_atoms: []})
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return <div className={clsx(className, styles.root)}>
    <NewPeriodicTable
      aggregations={statistics.atoms}
      metric={metric}
      exclusive={exclusive}
      values={elements}
      onChanged={handleElementsChanged}
      onExclusiveChanged={handleExclusiveChanged}
    />
    <Grid container spacing={2} className={styles.grid}>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.chemical_formula_hill"
          label="formula"
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.chemical_formula_anonymous"
          label="formula anonymous"
        />
      </Grid>
      <Grid item xs={6}>
        <FilterText
          quantity="results.material.n_elements"
          label="number of species"
        />
      </Grid>
    </Grid>
  </div>
})
FilterElements.propTypes = {
  className: PropTypes.string
}

export default FilterElements
