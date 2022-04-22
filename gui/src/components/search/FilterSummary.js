
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
import { isNil, isPlainObject } from 'lodash'
import { FilterChip, FilterChipGroup } from './FilterChip'
import { useSearchContext } from './SearchContext'
import { filterAbbreviations } from './FilterRegistry'
import { useUnits } from '../../units'
import { DType } from '../../utils'

/**
 * Displays a summary for the given subset of filters. Each filter value is
 * displayed as a chip.
 */
const typesWithLabel = new Set([DType.Boolean, DType.Int, DType.Float])
const useStyles = makeStyles(theme => {
  const paddingVertical = theme.spacing(1)
  const paddingHorizontal = theme.spacing(1.5)
  return {
    root: {
      boxShadow: 'inset 0 0 8px 1px rgba(0,0,0, 0.06)',
      backgroundColor: theme.palette.background.default,
      width: '100%',
      paddingRight: paddingHorizontal,
      paddingLeft: paddingHorizontal,
      paddingTop: paddingVertical,
      paddingBottom: paddingVertical,
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
  const { filterData, useFiltersState, useFiltersLockedState } = useSearchContext()
  const [filters, setFilter] = useFiltersState(quantities)
  const filtersLocked = useFiltersLockedState(quantities)
  const theme = useTheme()
  const units = useUnits()
  const styles = useStyles({classes: classes, theme: theme})
  let chips = []

  // Creates a set of chips for a quantity
  const createChips = useCallback((name, label, filterValue, onDelete) => {
    const newChips = []
    const locked = filtersLocked[name]
    if (isNil(filterValue)) {
      return
    }
    // Is query has multiple elements, we display a chip for each. For
    // numerical values we also display the quantity name.
    const metaType = filterData[name].dtype
    const serializer = filterData[name].serializerPretty
    const isArray = filterValue instanceof Array
    const isSet = filterValue instanceof Set
    const isObj = isPlainObject(filterValue)
    if (isArray || isSet) {
      filterValue.forEach((value, index) => {
        const displayValue = serializer(value, units)
        const item = <FilterChip
          key={`${name}${index}`}
          locked={locked}
          quantity={name}
          label={typesWithLabel.has(metaType) ? `${label} = ${displayValue}` : displayValue}
          onDelete={() => {
            let newValue
            if (isSet) {
              newValue = new Set(filterValue)
              newValue.delete(value)
            } else if (isArray) {
              newValue = [...filterValue]
              newValue.splice(index, 1)
            }
            onDelete(newValue)
          }}
        />
        newChips.push(item)
      })
    } else if (isObj) {
      // Range queries
      let lte = serializer(filterValue.lte, units)
      let gte = serializer(filterValue.gte, units)
      let lt = serializer(filterValue.lt, units)
      let gt = serializer(filterValue.gt, units)
      if (!isNil(gte) || !isNil(gt) || !isNil(lte) || !isNil(lt)) {
        let content
        if ((!isNil(gte) || !isNil(gt)) && (isNil(lte) && isNil(lt))) {
          content = `${label}${!isNil(gte) ? ` >= ${gte}` : ''}${!isNil(gt) ? ` > ${gt}` : ''}`
        } else if ((!isNil(lte) || !isNil(lt)) && (isNil(gte) && isNil(gt))) {
          content = `${label}${!isNil(lte) ? ` <= ${lte}` : ''}${!isNil(lt) ? ` < ${lt}` : ''}`
        } else {
          content = `${!isNil(gte) ? `${gte} <= ` : ''}${!isNil(gt) ? `${gt} < ` : ''}${label}${!isNil(lte) ? ` <= ${lte}` : ''}${!isNil(lt) ? ` < ${lt}` : ''}`
        }
        const item = <FilterChip
          key={name}
          locked={locked}
          quantity={name}
          label={content}
          onDelete={() => {
            onDelete(undefined)
          }}
        />
        newChips.push(item)
      }
    } else {
      const item = <FilterChip
        key={name}
        locked={locked}
        quantity={name}
        label={`${label}=${serializer(filterValue)}`}
        onDelete={() => {
          onDelete(undefined)
        }}
      />
      newChips.push(item)
    }
    return newChips
  }, [filterData, filtersLocked, units])

  if (quantities) {
    for (let quantity of quantities) {
      const filterValue = filters[quantity]
      if (isNil(filterValue)) {
        continue
      }
      // Nested filters
      const isSection = filterData[quantity].section
      let newChips
      if (isSection) {
        newChips = []
        for (let [key, value] of Object.entries(filterValue)) {
          const onDelete = (newValue) => {
            let newSection = {...filterValue}
            if (newValue === undefined) {
              delete newSection[key]
            } else {
              newSection[key] = newValue
            }
            setFilter([quantity, newSection])
          }
          newChips = newChips.concat(createChips(`${quantity}.${key}`, key, value, onDelete))
        }
      // Regular non-nested filters
      } else {
        const onDelete = (newValue) => setFilter([quantity, newValue])
        const label = filterAbbreviations[quantity]
        newChips = createChips(quantity, label, filterValue, onDelete)
      }
      if (newChips.length > 0) {
        const group = <FilterChipGroup
          key={quantity}
          quantity={quantity}
        >{newChips}
        </FilterChipGroup>
        chips.push(group)
      }
    }
  }

  return chips.length !== 0 && <div className={clsx(className, styles.root)}>
    {chips}
  </div>
})

FilterSummary.propTypes = {
  quantities: PropTypes.object, // Set of searchQuantities for which the filters are displayed
  className: PropTypes.string,
  classes: PropTypes.object
}

export default FilterSummary
