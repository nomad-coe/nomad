
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
import React, { useMemo } from 'react'
import { selector, useRecoilState } from 'recoil'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Chip } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { filterFamily } from './FilterContext'

/**
 * Displays an interactable summary for a given subset of filters
 * (=searchQuantities).
 */
const useStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
  return {
    root: {
      boxShadow: 'inset 0 0 8px 1px rgba(0,0,0, 0.075)',
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
  filters,
  id,
  className,
  classes
}) => {
  // We dynamically create a Recoil.js selector that is subscribed to the
  // filters specified in the input. This way only the specified filters will
  // cause a render.
  const filterState = useMemo(() => {
    return selector({
      key: id,
      get: ({get}) => {
        const query = {}
        for (let key of filters) {
          const filter = get(filterFamily(key))
          query[key] = filter
        }
        return query
      },
      set: ({set}, [key, value]) => {
        set(filterFamily(key), value)
      }
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  const [query, setQuery] = useRecoilState(filterState)
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const chips = []
  for (let filterName of filters) {
    const filterValue = query[filterName]
    const filterAbbr = filterName.split('.').pop()
    if (!filterValue) {
      continue
    }
    // Is query has multiple elements, we display a chip for each without
    // showing the quantity name
    const isArray = filterValue instanceof Array
    const isSet = filterValue instanceof Set
    const isObj = typeof filterValue === 'object' && filterValue !== null
    if (isArray || isSet) {
      filterValue.forEach((value, index) => {
        const item = <div
          key={chips.length}
          className={styles.chip}
        >
          <Chip
            label={value}
            onDelete={() => {
              if (isSet) {
                filterValue.delete(value)
                setQuery([filterName, filterValue])
              } else if (isArray) {
                filterValue.splice(index, 1)
                setQuery([filterName, filterValue])
              }
            }}
            color="primary"
          />
        </div>
        chips.push(item)
      })
    // Is query is an object, we display a customized view based on its contents.
    } else if (isObj) {
      const lte = filterValue.lte
      const gte = filterValue.gte
      let label = `${gte ? `${gte}<=` : ''}${filterAbbr}${lte ? `<=${lte}` : ''}`
      const item = <div
        key={chips.length}
        className={styles.chip}
      >
        <Chip
          label={label}
          onDelete={() => {
            setQuery([filterName, undefined])
          }}
          color="primary"
        />
      </div>
      chips.push(item)
    // If query is scalar-like, we show the quantity name together with the
    // value
    } else {
      const item = <div
        key={chips.length}
        className={styles.chip}
      >
        <Chip
          label={`${filterAbbr}=${filterValue}`}
          onDelete={() => {
            setQuery([filterName, undefined])
          }}
          color="primary"
        />
      </div>
      chips.push(item)
    }
  }
  return chips.length !== 0 && <div className={clsx(className, styles.root)}>
    {chips}
  </div>
})

FilterSummary.propTypes = {
  filters: PropTypes.arrayOf(PropTypes.string).isRequired, // Set of searchQuantities for which the filters are displayed
  id: PropTypes.string.isRequired, // Unique id for this summary. Needed in order to dynamically construct a Recoil.js Selector
  className: PropTypes.string,
  classes: PropTypes.object
}

export default FilterSummary
