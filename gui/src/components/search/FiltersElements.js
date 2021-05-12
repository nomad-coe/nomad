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
import React, { useState, useContext } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import PeriodicTable from './PeriodicTable'
import { searchContext } from './SearchContext'

const useFiltersElementStyles = makeStyles(theme => ({
  root: {}
}))

/**
 * Displays the filter options for chemical elements.
 */
const FiltersElements = React.memo(({
  className
}) => {
  const styles = useFiltersElementStyles()
  const [exclusive, setExclusive] = useState(false)
  const {response: {statistics, metric}, query, setQuery, setStatistics} = useContext(searchContext)

  const handleExclusiveChanged = () => {
    const newExclusive = !exclusive
    if (newExclusive) {
      setQuery({only_atoms: query.atoms, atoms: []})
    } else {
      setQuery({atoms: query.only_atoms, only_atoms: []})
    }
    setExclusive(newExclusive)
  }

  const handleAtomsChanged = atoms => {
    exclusive && setExclusive(false)
    setQuery({atoms: atoms, only_atoms: []})
  }

  return <div className={clsx(className, styles.root)}>
    <PeriodicTable
      aggregations={statistics.atoms}
      metric={metric}
      exclusive={exclusive}
      values={[...(query.atoms || []), ...(query.only_atoms || [])]}
      onChanged={handleAtomsChanged}
      onExclusiveChanged={handleExclusiveChanged}
    />
  </div>
})
FiltersElements.propTypes = {
  className: PropTypes.string
}

export default FiltersElements
