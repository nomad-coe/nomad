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
import { Box } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useRecoilValue } from 'recoil'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import InputItem, { inputItemHeight } from './InputItem'
import InputUnavailable from './InputUnavailable'
import Placeholder from '../../visualization/Placeholder'
import { useSearchContext } from '../SearchContext'
import { isNil, isNumber } from 'lodash'
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
    boxSizing: 'border-box',
    marginBottom: theme.spacing(-0.5)
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
const InputTerms = React.memo(({
  searchQuantity,
  title,
  description,
  visible,
  nColumns,
  scale,
  increment,
  showInput,
  showHeader,
  showStatistics,
  showSuggestions,
  options,
  sortStatic,
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
  const [scaleState, setScaleState] = useState(scale || filterData[searchQuantity]?.scale || 'linear')
  showStatistics = isStatisticsEnabled && showStatistics
  const nOptions = isNumber(options) ? options : undefined

  // See if the filter has a fixed amount of options. These may have been
  // explicitly provided or defined in the metainfo. If options is given as a
  // number, any fixed options are ignored and the data is retrieved through an
  // aggregation.
  const fixedOptions = useMemo(() => isNil(nOptions)
      ? options || filterData[searchQuantity]?.options
      : undefined
  , [nOptions, filterData, searchQuantity, options])

  const nFixedOptions = fixedOptions && Object.keys(fixedOptions).length
  const minSize = nOptions === 0
    ? 0
    : nOptions || nFixedOptions || filterData[searchQuantity]?.aggs?.terms?.size || 5
  const placeholder = filterData[searchQuantity]?.placeholder || "Type here"
  const [requestedAggSize, setRequestedAggSize] = useState(minSize)
  const incr = useState(increment || minSize)[0]
  const [loading, setLoading] = useState(false)

  // If a fixed list of options is used, we must restrict the aggregation return
  // values with 'include'. Otherwise the returned results may contain other
  // values.
  const aggConfig = useMemo(() => {
    const config = {type: 'terms', size: minSize}
    if (fixedOptions) config.include = Object.keys(fixedOptions)
    return config
  }, [minSize, fixedOptions])

  const agg = useAgg(searchQuantity, visible && !(nOptions === 0) && !(!showStatistics && fixedOptions), 'scroll', aggConfig)
  const aggCall = useAggCall(searchQuantity, 'scroll')
  const receivedAggSize = agg?.data?.length
  const [filter, setFilter] = useFilterState(searchQuantity)
  const unavailable = (nOptions === 0) ? false : !(agg?.data && agg.data.length > 0)
  const disabled = unavailable

  // Form the final list of options. If no fixed options are available, the
  // options are gathered from the aggregation.
  const finalOptions = useMemo(() => {
    if (fixedOptions) return fixedOptions
    if (!agg?.data) return {}

    const maxSize = Math.min(requestedAggSize, agg.data.length)
    return agg.data.slice(0, maxSize).reduce((opt, { value }) => {
      opt[value] = { label: value }
      return opt
    }, {})
  }, [fixedOptions, agg, requestedAggSize])

  // Modify the checkboxes according to changing filters, changing aggregation
  // results or change in the available options.
  useEffect(() => {
    let options = Object.entries(finalOptions).reduce((opt, [key, value]) => {
      const selected = filter?.has(key) || false
      opt[key] = {
        checked: selected,
        label: value.label,
        disabled: isStatisticsEnabled && showStatistics && !selected
      }
      return opt
    }, {})

    if (agg?.data) {
      // Update counts and disable if not selected and count is 0
      agg.data?.forEach(({ value, nested_count }) => {
        const selected = filter?.has(value) || false
        if (options[value]) {
          options[value].nested_count = nested_count
          options[value].disabled = selected ? false : nested_count === 0
        }
      })

      // Sort by count if using fixed options and sorting is enabled
      if (fixedOptions && sortStatic) {
        options = Object.fromEntries(
          Object.entries(options).sort(([, a], [, b]) => {
            const bCount = b.nested_count || 0
            const aCount = a.nested_count || 0
            return bCount - aCount
          })
        )
      }
    }

    setVisibleOptions(options)
  }, [agg?.data, filter, finalOptions, fixedOptions, isStatisticsEnabled, showStatistics, sortStatic])

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
    return showInput
      ? <InputTooltip unavailable={unavailable}>
        <div className={styles.container}>
          <InputTextQuantity
            className={styles.textField}
            quantity={searchQuantity}
            disabled={disabled}
            disableSuggestions={!showSuggestions}
            placeholder={placeholder}
            fullWidth
          />
        </div>
      </InputTooltip>
      : null
  }, [showInput, unavailable, styles, searchQuantity, disabled, showSuggestions, placeholder])

  // Create the options component
  const optionsComponent = useMemo(() => {
    if ((nOptions === 0)) {
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

    const xs = 12 / nColumns
    const nRows = Math.ceil(nItems * xs / 12)
    const actionsHeight = 34
    let reservedHeight
    if (aggCollapse === 'on') {
      const itemHeight = nRows * inputItemHeight
      const actionHeight = (showMore || showLess) ? actionsHeight : 0
      reservedHeight = nItems > 0 ? (itemHeight + actionHeight) : undefined
    } else if (aggCollapse === 'off') {
      const itemHeight = Math.max(minSize * xs / 12, nRows) * inputItemHeight
      const allLoaded = !isNil(nFixedOptions) && nOptions >= nFixedOptions
      const actionHeight = (fixedOptions || allLoaded) ? 0 : actionsHeight
      reservedHeight = itemHeight + actionHeight
    }

    const max = agg ? Math.max(...agg.data.map(option => option.nested_count)) : 0
    const items = visibleOptions && <div
      className={styles.grid}
      style={{gridTemplateRows: `repeat(${nRows}, 1fr)`}}
    >
      {Object.entries(visibleOptions).map(([key, value]) => (
        <InputItem
          key={key}
          value={key}
          tooltip={value.description}
          label={value.label}
          selected={value.checked}
          disabled={value.disabled}
          disableStatistics={!showStatistics}
          onChange={handleChange}
          variant="checkbox"
          max={max}
          count={value.nested_count}
          scale={scaleState}
        />
      ))}
    </div>

    const noMore = agg?.exhausted && receivedAggSize === requestedAggSize
    let aggComp

    // If statistics are disabled and there is a fixed number of options, we
    // show the options immediately.
    if ((!showStatistics || !sortStatic) && fixedOptions) {
      aggComp = items
    // No fixed options or aggregation data is available
    } else if (receivedAggSize === 0 && !fixedOptions) {
      aggComp = <InputUnavailable/>
    // Show placeholder if aggregation is underway
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
    showStatistics,
    agg,
    finalOptions,
    minSize,
    fixedOptions,
    receivedAggSize,
    requestedAggSize,
    incr,
    nColumns,
    aggCollapse,
    visibleOptions,
    styles,
    aggIndicator,
    handleChange,
    scaleState,
    handleShowMore,
    handleShowLess,
    loading,
    nFixedOptions,
    testID,
    nOptions,
    sortStatic
  ])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    {showHeader
      ? <InputHeader
        quantity={searchQuantity}
        label={title}
        description={description}
        scale={scaleState}
        onChangeScale={setScaleState}
        disableStatistics={!showStatistics}
      />
      : <Box marginBottom={0.5}/>
    }
    {searchComponent}
    {optionsComponent}
  </div>
})

InputTerms.propTypes = {
  searchQuantity: PropTypes.string.isRequired,
  title: PropTypes.string,
  description: PropTypes.string,
  visible: PropTypes.bool,
  nColumns: PropTypes.number, // Number of columns to use
  scale: PropTypes.string, // The statistics scaling
  options: PropTypes.oneOfType([PropTypes.number, PropTypes.object]), // Controls what options to show
  sortStatic: PropTypes.bool, // Whether to sort statically defined options by occurrence
  increment: PropTypes.number, // The amount of new items to load on 'show more'
  showInput: PropTypes.bool, // Whether to show the search input field
  showHeader: PropTypes.bool, // Whether to show the header
  showStatistics: PropTypes.bool, // Whether to disable statistics
  showSuggestions: PropTypes.bool, // Whether to disable the text field suggestions
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

InputTerms.defaultProps = {
  nColumns: 1,
  showHeader: true,
  showInput: true,
  showStatistics: true,
  sortStatic: true,
  'data-testid': 'input-terms'
}

export default InputTerms
