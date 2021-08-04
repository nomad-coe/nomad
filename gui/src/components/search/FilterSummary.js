
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
import React, { useCallback } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import FilterChip from './FilterChip'
import { useFiltersState } from './FilterContext'
import { formatNumber, getIsNumeric } from '../../utils'
import { Quantity, useUnits } from '../../units'

/**
 * Displays an interactable summary for a given subset of filters
 * (=searchQuantities).
 */
const useStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
  return {
    root: {
      boxShadow: 'inset 0 0 8px 1px rgba(0,0,0, 0.06)',
      backgroundColor: theme.palette.background.default,
      width: '100%',
      padding: padding,
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'row',
      flexWrap: 'wrap'
    },
    chip: {
      padding: theme.spacing(0.5)
    }
  }
})
const FilterSummary = React.memo(({
  quantities,
  className,
  classes
}) => {
  const [filters, setFilter] = useFiltersState(quantities)
  const theme = useTheme()
  const units = useUnits()
  const styles = useStyles({classes: classes, theme: theme})
  const chips = []

  // Formats values so that they are more readable.
  const format = useCallback((value) => {
    if (!isNil(value)) {
      if (value instanceof Quantity) {
        value = `${formatNumber(value.toSystem(units))} ${value.unit.label(units)}`
      } else {
        value = formatNumber(value)
      }
    }
    return value
  }, [units])

  for (let quantity of quantities) {
    const filterValue = filters[quantity]
    const filterAbbr = quantity.split('.').pop()
    if (!filterValue) {
      continue
    }
    const isNumeric = getIsNumeric(quantity)
    // Is query has multiple elements, we display a chip for each. For numerical
    // values we also display the quantity name.
    const isArray = filterValue instanceof Array
    const isSet = filterValue instanceof Set
    const isObj = typeof filterValue === 'object' && filterValue !== null
    if (isArray || isSet) {
      filterValue.forEach((value, index) => {
        const displayValue = format(value)
        const item = <FilterChip
          key={chips.length}
          label={isNumeric ? `${filterAbbr} = ${displayValue}` : displayValue}
          onDelete={() => {
            if (isSet) {
              const newSet = new Set(filterValue)
              newSet.delete(value)
              setFilter([quantity, newSet])
            } else if (isArray) {
              const newArray = [...filterValue]
              newArray.splice(index, 1)
              setFilter([quantity, newArray])
            }
          }}
        />
        chips.push(item)
      })
    // Is query is an object, it is assumed to represent a range filter.
    } else if (isObj) {
      let lte = format(filterValue.lte)
      let gte = format(filterValue.gte)
      let lt = format(filterValue.lt)
      let gt = format(filterValue.gt)
      let label
      if ((!isNil(gte) || !isNil(gt)) && (isNil(lte) && isNil(lt))) {
        label = `${filterAbbr}${!isNil(gte) ? ` >= ${gte}` : ''}${!isNil(gt) ? ` > ${gt}` : ''}`
      } else if ((!isNil(lte) || !isNil(lt)) && (isNil(gte) && isNil(gt))) {
        label = `${filterAbbr}${!isNil(lte) ? ` <= ${lte}` : ''}${!isNil(lt) ? ` < ${lt}` : ''}`
      } else {
        label = `${!isNil(gte) ? `${gte} <= ` : ''}${!isNil(gt) ? `${gt} < ` : ''}${filterAbbr}${!isNil(lte) ? ` <= ${lte}` : ''}${!isNil(lt) ? ` < ${lt}` : ''}`
      }
      const item = <FilterChip
        key={chips.length}
        label={label}
        onDelete={() => {
          setFilter([quantity, undefined])
        }}
      />
      chips.push(item)
    // If query is scalar-like, we show the quantity name together with the
    // value
    } else {
      const item = <FilterChip
        key={chips.length}
        label={`${filterAbbr}=${filterValue}`}
        onDelete={() => {
          setFilter([quantity, undefined])
        }}
      />
      chips.push(item)
    }
  }
  return chips.length !== 0 && <div className={clsx(className, styles.root)}>
    {chips}
  </div>
})

FilterSummary.propTypes = {
  quantities: PropTypes.arrayOf(PropTypes.string).isRequired, // Set of searchQuantities for which the filters are displayed
  className: PropTypes.string,
  classes: PropTypes.object
}

export default FilterSummary
