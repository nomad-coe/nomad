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
import React, { useCallback, useEffect, useState, useRef } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Grid,
  FormControlLabel,
  Checkbox
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
import { useFilterState, useAgg } from './FilterContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  }
}))
const FilterCheckboxes = React.memo(({
  label,
  quantity,
  description,
  options,
  visible,
  xs,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [visibleOptions, setVisibleOptions] = useState()
  const availableOptions = useAgg(quantity, 'terms', true, visible)
  const [filter, setFilter] = useFilterState(quantity)
  const firstFetch = useRef(true)

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const title = name

  // Save the available options when retrieved for the first time (without any
  // filters)
  useEffect(() => {
    if (availableOptions && firstFetch.current) {
      const opt = {}
      for (let option of availableOptions) {
        if (option.count > 0) {
          if (!options || options[option.value]) {
            opt[option.value] = {
              checked: false,
              disabled: false,
              label: options && options[option.value].label
            }
          }
        }
      }
      setVisibleOptions(opt)
      firstFetch.current = false
    }
  }, [availableOptions, options])

  // React to changing filters
  useEffect(() => {
    setVisibleOptions(old => {
      const newOptions = {}
      if (old === undefined) {
        return old
      }
      for (let key of Object.keys(old)) {
        newOptions[key] = old[key]
        newOptions[key].checked = filter ? filter.has(key) : false
      }
      return newOptions
    })
  }, [filter])

  // React to changing options. Unselected ones will be disabled/enabled
  // accordingly.
  useEffect(() => {
    if (availableOptions) {
      setVisibleOptions(old => {
        const newOptions = {}
        for (let key of Object.keys(old)) {
          newOptions[key] = old[key]
          if (!newOptions[key].checked) {
            newOptions[key].disabled = availableOptions[key]?.count > 0
          }
        }
        return newOptions
      })
    }
  }, [availableOptions])

  const handleChange = useCallback((event) => {
    const newOptions = {...visibleOptions}
    newOptions[event.target.name].checked = event.target.checked
    const checked = Object.entries(newOptions)
      .filter(([key, value]) => value.checked)
      .map(([key, value]) => key)
    setFilter(new Set(checked))
  }, [setFilter, visibleOptions])

  const checkboxes = visibleOptions && Object.entries(visibleOptions).map(([key, value]) => {
    return <Grid item xs={xs} key={key}>
      <FormControlLabel
        control={<Checkbox checked={value.checked} onChange={handleChange} name={key}/>}
        label={value.label || key}
        disabled={value.disabled}
      />
    </Grid>
  })

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <FilterLabel label={title} description={desc}/>
    <Grid container spacing={0}>
      {checkboxes}
    </Grid>
  </div>
})

FilterCheckboxes.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  // Optional information about the options. Can also be used to enable/disable
  // options.
  options: PropTypes.object,
  visible: PropTypes.bool,
  xs: PropTypes.number,
  description: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

FilterCheckboxes.defaultProps = {
  xs: 12
}

export default FilterCheckboxes
