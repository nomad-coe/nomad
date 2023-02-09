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
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useRecoilValue } from 'recoil'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import InputItem, { inputItemHeight } from './InputItem'
import InputUnavailable from './InputUnavailable'
import Placeholder from '../../visualization/Placeholder'
import { formatLabel } from '../../../utils'
import { useSearchContext } from '../SearchContext'
import { isNil } from 'lodash'
import Pagination from '../../visualization/Pagination'
import { guiState } from '../../GUIMenu'
import { InputTextQuantity } from './InputText'

/**
 * Generic input component that can be configured to fit several use cases. The
 * most typical configufations are:
 * - Smallish, fixed number of options: do not show search field, use
 *   metainfo/explicit options to always show the options.
 * - Large amount of options which are not unique to entries: show search field
 *   and a subset of options sorted by their occurrence.
 * - Large amount of options which are unique for each entry: only show search field.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  container: {
    width: '100%'
  },
  grid: {
    width: '100%',
    display: 'grid',
    gridAutoFlow: 'column'
  },
  progress: {
    marginLeft: theme.spacing(0.5)
  },
  textField: {
    marginBottom: theme.spacing(1)
  }
}))
const InputField = React.memo(({
  quantity,
  label,
  description,
  visible,
  xs,
  initialScale,
  initialSize,
  increment,
  disableStatistics,
  disableSearch,
  disableOptions,
  disableSuggestions,
  formatLabels,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const {
    filterData,
    useAgg,
    useAggCall,
    useFilterState,
    useIsStatisticsEnabled
  } = useSearchContext()
  const isStatisticsEnabled = useIsStatisticsEnabled()
  const styles = useStyles({classes: classes, theme: theme})
  const [visibleOptions, setVisibleOptions] = useState()
  const aggIndicator = useRecoilValue(guiState('aggIndicator'))
  const aggCollapse = useRecoilValue(guiState('aggCollapse'))
  const [scale, setScale] = useState(initialScale || filterData[quantity]?.scale || 'linear')
  disableStatistics = isNil(disableStatistics) ? !isStatisticsEnabled : disableStatistics

  // See if the filter has a fixed amount of options. These may have been
  // explicitly provided or defined in the metainfo. If you explicitly specify
  // an initialSize, any fixed options are ignored and the data is retrieved
  // through an aggregation (this is done because the top aggregations may not
  // match the list of explicit options).
  const fixedOptions = useMemo(() => {
    return isNil(initialSize)
      ? filterData[quantity]?.options
      : undefined
  }, [initialSize, filterData, quantity])
  const nFixedOptions = fixedOptions && Object.keys(fixedOptions).length
  const minSize = disableOptions ? 0 : initialSize || nFixedOptions || filterData[quantity]?.aggs?.terms?.size
  const placeholder = filterData[quantity]?.placeholder || "Type here"
  const [requestedAggSize, setRequestedAggSize] = useState(minSize)
  const incr = useState(increment || minSize)[0]
  const [loading, setLoading] = useState(false)
  const aggConfig = useMemo(() => {
    const config = {type: 'terms', size: minSize}
    // If a fixed list of options is used, we must restrict the aggregation
    // return values with 'include'. Otherwise the returned results may not
    // contain the correct values.
    const options = filterData[quantity]?.options
    if (options) config.include = Object.keys(options)
    return config
  }, [minSize, filterData, quantity])
  const agg = useAgg(quantity, visible && !disableOptions && !(disableStatistics && fixedOptions), 'scroll', aggConfig)
  const aggCall = useAggCall(quantity, 'scroll')
  const receivedAggSize = agg?.data?.length
  const [filter, setFilter] = useFilterState(quantity)
  const unavailable = disableOptions ? false : !(agg?.data && agg.data.length > 0)
  const disabled = unavailable

  // Form the final list of options. If no fixed options are available, the
  // options are gathered from the aggregation.
  const finalOptions = useMemo(() => {
    if (fixedOptions) {
      return fixedOptions
    }
    if (agg?.data) {
      const opt = {}
      const maxSize = Math.min(requestedAggSize, agg.data.length)
      for (let i = 0; i < maxSize; ++i) {
        const value = agg.data[i]
        opt[value.value] = {label: value.value}
      }
      return opt
    }
    return {}
  }, [fixedOptions, agg, requestedAggSize])

  // Modify the checkboxes according to changing filters, changing aggregation
  // results or change in the available options.
  useEffect(() => {
    const opt = {}
    for (const [key, value] of Object.entries(finalOptions)) {
      opt[key] = {
        checked: filter ? filter.has(key) : false,
        label: formatLabels ? formatLabel(value.label) : value.label,
        disabled: !disableStatistics
      }
    }
    if (agg?.data) {
      for (const value of agg.data) {
        const key = value.value
        const selected = filter ? filter.has(key) : false
        const oldState = opt[key]
        const disabled = selected ? false : value.count === 0
        if (oldState) {
          oldState.count = value.count
          oldState.disabled = disabled
        }
      }
    }
    setVisibleOptions(opt)
  }, [agg, filter, finalOptions, disableStatistics, formatLabels])

  // Show more values
  const handleShowMore = useCallback(() => {
    setLoading(true)
    let newSize = requestedAggSize + incr
    aggCall({type: 'terms', size: newSize}, (response) => {
      if (response?.data?.length === requestedAggSize) {
        newSize = requestedAggSize
      }
      setLoading(false)
      setRequestedAggSize(newSize)
    })
  }, [aggCall, incr, requestedAggSize])

  // Show less values
  const handleShowLess = useCallback(() => {
    setRequestedAggSize(old => {
      const newSize = Math.max(old - incr, minSize)
      aggCall({type: 'terms', size: newSize}, () => {})
      return newSize
    })
  }, [aggCall, incr, minSize])

  // Handle changes in the selected values
  const handleChange = useCallback((event, key, selected) => {
    const newOptions = {...visibleOptions}
    newOptions[key].checked = selected
    const checked = Object.entries(newOptions)
      .filter(([key, value]) => value.checked)
      .map(([key, value]) => key)
    setFilter(new Set(checked))
  }, [setFilter, visibleOptions])

  // Create the search component
  const searchComponent = useMemo(() => {
    return disableSearch
      ? null
      : <InputTooltip unavailable={unavailable}>
        <div className={styles.container}>
          <InputTextQuantity
            className={styles.textField}
            quantity={quantity}
            disabled={disabled}
            disableSuggestions={disableSuggestions}
            placeholder={placeholder}
            fullWidth
          />
        </div>
      </InputTooltip>
  }, [disableSearch, unavailable, styles, quantity, disabled, disableSuggestions, placeholder])

  // Create the options component
  const optionsComponent = useMemo(() => {
    if (disableOptions) {
      return
    }

    const nItems = agg ? Object.keys(finalOptions).length : minSize
    const didNotReceiveMore = isNil(receivedAggSize) ? false : receivedAggSize < requestedAggSize
    const noMoreAvailable = isNil(nFixedOptions) ? false : requestedAggSize >= nFixedOptions
    const hide = didNotReceiveMore || noMoreAvailable
    const showMore = fixedOptions
      ? false
      : !hide
    const showLess = fixedOptions
      ? false
      : (requestedAggSize - incr >= minSize)

    const nRows = Math.ceil(nItems * xs / 12)
    const actionsHeight = 34
    let reservedHeight
    if (aggCollapse === 'on') {
      const itemHeight = nRows * inputItemHeight
      const actionHeight = (showMore || showLess) ? actionsHeight : 0
      reservedHeight = nItems > 0 ? (itemHeight + actionHeight) : undefined
    } else if (aggCollapse === 'off') {
      const itemHeight = Math.max(minSize * xs / 12, nRows) * inputItemHeight
      const allLoaded = !isNil(nFixedOptions) && initialSize >= nFixedOptions
      const actionHeight = (fixedOptions || allLoaded) ? 0 : actionsHeight
      reservedHeight = itemHeight + actionHeight
    }

    const max = agg ? Math.max(...agg.data.map(option => option.count)) : 0
    const items = visibleOptions && <div
      className={styles.grid}
      style={{gridTemplateRows: `repeat(${nRows}, 1fr)`}}
    >
      {Object.entries(visibleOptions).map(([key, value]) => (
        <InputItem
          key={key}
          value={key}
          label={value.label}
          selected={value.checked}
          disabled={value.disabled}
          disableStatistics={disableStatistics}
          onChange={handleChange}
          variant="checkbox"
          max={max}
          count={value.count}
          scale={scale}
        />
      ))}
    </div>

    const noMore = agg?.exhausted && receivedAggSize === requestedAggSize
    let aggComp
    if (fixedOptions) {
      aggComp = items
    } else if (receivedAggSize === 0) {
      aggComp = <InputUnavailable/>
    } else if (!agg && aggIndicator === 'on') {
      aggComp = <Placeholder
        variant="rect"
        data-testid={`${testID}-placeholder`}
        margin={0}
      />
    } else {
      aggComp = <>
        {items}
        <Pagination
          showMore={showMore}
          showLess={showLess}
          disableMore={noMore}
          loadingMore={loading}
          marginTop={0.5}
          onLess={handleShowLess}
          onMore={handleShowMore}
        />
      </>
    }
    return <div className={styles.container} style={{height: reservedHeight}}>
      {aggComp}
    </div>
  }, [
    disableOptions,
    disableStatistics,
    agg,
    finalOptions,
    minSize,
    fixedOptions,
    receivedAggSize,
    requestedAggSize,
    incr,
    xs,
    aggCollapse,
    visibleOptions,
    styles,
    aggIndicator,
    initialSize,
    handleChange,
    scale,
    handleShowMore,
    handleShowLess,
    loading,
    nFixedOptions,
    testID
  ]
  )

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <InputHeader
      quantity={quantity}
      label={label}
      description={description}
      scale={scale}
      onChangeScale={setScale}
      disableStatistics={disableStatistics}
    />
    {searchComponent}
    {optionsComponent}
  </div>
})

InputField.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  visible: PropTypes.bool,
  xs: PropTypes.number,
  initialScale: PropTypes.number, // The initial statistics scaling
  initialSize: PropTypes.number, // The initial maximum number of items to load
  increment: PropTypes.number, // The amount of new items to load on 'show more'
  disableStatistics: PropTypes.bool, // Whether to disable statistics
  disableSearch: PropTypes.bool, // Whether to show the search field
  disableOptions: PropTypes.bool, // Whether to show the options gathered through aggregations
  disableSuggestions: PropTypes.bool, // Whether to disable the text field suggestions
  formatLabels: PropTypes.bool, // Whether to reformat the options labels
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

InputField.defaultProps = {
  xs: 12
}

export default InputField
