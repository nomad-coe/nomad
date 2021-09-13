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
import React, { useCallback, useEffect, useState, useMemo } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Grid,
  FormControlLabel,
  Checkbox
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../../searchQuantities'
import InputLabel from './InputLabel'
import { useFilterState, useAgg, useInitialAgg } from '../FilterContext'
import { isArray } from 'lodash'

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
const InputCheckboxes = React.memo(({
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
  const availableOptions = useAgg(quantity, true, visible)
  const initialAgg = useInitialAgg(quantity)
  const [filter, setFilter] = useFilterState(quantity)
  const finalOptions = useMemo(() => {
    // If explicit options provided, use them
    if (options) {
      return options
    }
    // If the metainfo has enumerable options, use them
    const quantities = searchQuantities?.[quantity]?.type?.type_data
    if (isArray(quantities) && quantities.length > 0) {
      const opt = {}
      for (const name of quantities) {
        opt[name] = {label: name}
      }
      return opt
    }
    // As a last resort, use the initial aggregation results as options.
    if (initialAgg) {
      const opt = {}
      for (const value of initialAgg) {
        opt[value.value] = {label: value.value}
      }
      return opt
    }
    return {}
  }, [options, quantity, initialAgg])

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const title = name

  // Modify the checkboxes according to changing filters, changing aggregation
  // results or change in the available options.
  useEffect(() => {
    const opt = {}
    for (let [key, value] of Object.entries(finalOptions)) {
      opt[key] = {
        checked: filter ? filter.has(key) : false,
        label: value.label,
        disabled: true
      }
    }
    if (availableOptions) {
      for (let value of availableOptions) {
        const key = value.value
        const selected = filter ? filter.has(key) : false
        const oldState = opt[key]
        const disabled = selected ? false : value.count === 0
        if (oldState) {
          oldState.disabled = disabled
        }
      }
    }
    setVisibleOptions(opt)
  }, [availableOptions, filter, finalOptions])

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
    <InputLabel label={title} description={desc}/>
    <Grid container spacing={0}>
      {checkboxes}
    </Grid>
  </div>
})

InputCheckboxes.propTypes = {
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

InputCheckboxes.defaultProps = {
  xs: 12
}

export default InputCheckboxes
