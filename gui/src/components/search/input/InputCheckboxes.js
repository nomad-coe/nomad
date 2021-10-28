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
import { Grid } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../../searchQuantities'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import InputItem from './InputItem'
import {
  useFilterState,
  useAgg,
  useInitialAgg,
  useFilterLocked,
  filterData
} from '../SearchContext'
import { isArray } from 'lodash'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  gridItem: {
    position: 'relative'
  }
}))
const InputCheckboxes = React.memo(({
  quantity,
  label,
  description,
  visible,
  xs,
  initialScale,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [visibleOptions, setVisibleOptions] = useState()
  const [scale, setScale] = useState(initialScale)
  const agg = useAgg(quantity, visible)
  const initialAgg = useInitialAgg(quantity)
  const [filter, setFilter] = useFilterState(quantity)
  const locked = useFilterLocked(quantity)
  const finalOptions = useMemo(() => {
    // If explicit options provided, use them
    const options = filterData[quantity].options
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
    // As a last resort, use the initial aggregation results as options. We
    // cannot use the agg results because the page may be loaded with additional
    // filters e.g. from the query string which will alter the available
    // options.
    if (initialAgg?.data) {
      const opt = {}
      for (const value of initialAgg.data) {
        opt[value.value] = {label: value.value}
      }
      return opt
    }
    return {}
  }, [quantity, initialAgg])

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
    if (agg?.data) {
      for (let value of agg.data) {
        const key = value.value
        const selected = filter ? filter.has(key) : false
        const oldState = opt[key]
        const disabled = locked || (selected ? false : value.count === 0)
        if (oldState) {
          oldState.count = value.count
          oldState.disabled = disabled
        }
      }
    }
    setVisibleOptions(opt)
  }, [agg, filter, finalOptions, locked])

  const handleChange = useCallback((key, value) => {
    const newOptions = {...visibleOptions}
    newOptions[key].checked = value
    const checked = Object.entries(newOptions)
      .filter(([key, value]) => value.checked)
      .map(([key, value]) => key)
    setFilter(new Set(checked))
  }, [setFilter, visibleOptions])

  const checkboxes = visibleOptions && Object.entries(visibleOptions).map(([key, value]) => {
    return <Grid item xs={xs} key={key} className={styles.gridItem}>
      <InputItem
        value={key}
        label={value.label}
        selected={value.checked}
        disabled={value.disabled}
        onChange={handleChange}
        variant="checkbox"
        total={agg?.total}
        count={value.count}
        scale={scale}
      />
    </Grid>
  })

  return <InputTooltip locked={locked} disabled={false}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <InputHeader
        quantity={quantity}
        label={label}
        description={description}
        scale={scale}
        onChangeScale={setScale}
        disableAggSize
      />
      <Grid container spacing={0}>
        {checkboxes}
      </Grid>
    </div>
  </InputTooltip>
})

InputCheckboxes.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  // Optional information about the options. Can also be used to enable/disable
  // options.
  visible: PropTypes.bool,
  xs: PropTypes.number,
  initialScale: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

InputCheckboxes.defaultProps = {
  xs: 12,
  initialScale: 1
}

export default InputCheckboxes
