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
import React, { useState, useMemo, useCallback, useEffect, useRef, useContext } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import InputHeader from './InputHeader'
import { inputSectionContext } from './InputSection'
import { Quantity } from '../../units/Quantity'
import { Unit } from '../../units/Unit'
import { useUnitContext } from '../../units/UnitContext'
import { DType, formatNumber } from '../../../utils'
import { getInterval } from '../../plotting/common'
import { useSearchContext } from '../SearchContext'
import PlotHistogram from '../../plotting/PlotHistogram'
import { isValid, getTime } from 'date-fns'
import { ActionCheckbox } from '../../Actions'
import { autorangeDescription } from '../widgets/WidgetHistogram'

/*
 * Component for displaying a numerical range as a slider/histogram together
 * with query input possibility.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%'
  },
  histogram: {
  }
}))
export const Range = React.memo(({
  xAxis,
  nSteps,
  visible,
  scale,
  nBins,
  disableHistogram,
  disableXTitle,
  autorange,
  showinput,
  aggId,
  className,
  classes,
  'data-testid': testID
}) => {
  const {filterData, useAgg, useFilterState, useIsStatisticsEnabled} = useSearchContext()
  const sectionContext = useContext(inputSectionContext)
  const repeats = sectionContext?.repeats
  const isStatisticsEnabled = useIsStatisticsEnabled()
  const styles = useStyles({classes})
  const [filter, setFilter] = useFilterState(xAxis.quantity)
  const [minLocal, setMinLocal] = useState()
  const [maxLocal, setMaxLocal] = useState()
  const [plotData, setPlotData] = useState({xAxis})
  const loading = useRef(false)
  const firstRender = useRef(true)
  const validRange = useRef()
  const rangeChanged = useRef(false)
  const minInputChanged = useRef(false)
  const maxInputChanged = useRef(false)
  const [range, setRange] = useState({gte: undefined, lte: undefined})
  const [minError, setMinError] = useState(false)
  const [maxError, setMaxError] = useState(false)
  const [minInclusive, setMinInclusive] = useState(true)
  const [maxInclusive, setMaxInclusive] = useState(true)
  const highlight = Boolean(filter)
  disableHistogram = isNil(disableHistogram) ? !isStatisticsEnabled : disableHistogram

  // Determine the description and units
  const def = filterData[xAxis.quantity]
  const unitStorage = useMemo(() => { return new Unit(def?.unit || 'dimensionless') }, [def])
  const discretization = xAxis.dtype === DType.Int ? 1 : undefined
  const isTime = xAxis.dtype === DType.Timestamp
  const firstLoad = useRef(true)

  // We need to set a valid initial input state: otherwise the component thinks
  // it is uncontrolled.
  const [minInput, setMinInput] = useState(isTime ? new Date() : '')
  const [maxInput, setMaxInput] = useState(isTime ? new Date() : '')

  // Function for converting filter values into storage units used by the API.
  const fromDisplayUnit = useCallback(filter => {
    return filter instanceof Quantity
      ? filter.to(unitStorage).value()
      : filter
  }, [unitStorage])

  // Function for converting filter values from storage unit to display unit
  const fromStorageUnit = useCallback(filter => {
    return (isTime)
      ? filter
      : filter instanceof Quantity
        ? filter.to(unitStorage)
        : new Quantity(filter, unitStorage)
  }, [unitStorage, isTime])

  // Function for converting filter values from display unit to storage unit
  const toStorageUnit = useCallback(filter => {
    return (isTime)
      ? filter
      : filter instanceof Quantity
        ? filter.to(unitStorage)
        : new Quantity(filter, xAxis.unit).to(unitStorage)
  }, [isTime, xAxis.unit, unitStorage])

  // Aggregation when the statistics are enabled: a histogram aggregation with
  // extended bounds based on the currently set filter range. Note: the config
  // should be memoed in order to prevent re-renders.
  const minRef = useRef(fromDisplayUnit(isNil(filter?.gte) ? filter?.gt : filter?.gte))
  const maxRef = useRef(fromDisplayUnit(isNil(filter?.lte) ? filter?.lt : filter?.lte))
  const aggHistogramConfig = useMemo(() => {
    const filterBounds = (filter) ? {
      min: fromDisplayUnit(isNil(filter.gte) ? filter.gt : filter.gte),
      max: fromDisplayUnit(isNil(filter.lte) ? filter.lt : filter.lte)
    } : undefined
    let exclude_from_search
    let extended_bounds
    // Extended bounds are always applied on the first load but with no
    // filtering.
    if (firstLoad.current) {
      extended_bounds = filterBounds
      exclude_from_search = true
    } else {
      // Is zoom is not disabled, the current filter defines the histogram
      // boundary.
      if (autorange) {
        extended_bounds = filterBounds
        exclude_from_search = false
      } else {
        // Is zoom is disabled and filter is outside the old boundary, the
        // filter defines the new boundary and filter is active.
        const outOfBounds = filter ? (filterBounds.min < minRef.current || filterBounds.max > maxRef.current) : false
        if (outOfBounds) {
          extended_bounds = filterBounds
          exclude_from_search = false
        } else {
          // Is zoom is disabled and filter is within the old boundary, and the
          // global boundary, the old boundary is kept.
          extended_bounds = filter
            ? (isNil(minRef.current) ? undefined : {min: minRef.current, max: maxRef.current})
            : undefined
          exclude_from_search = true
        }
      }
    }
    return (isTime || !discretization)
      ? {type: 'histogram', buckets: nBins, exclude_from_search, extended_bounds}
      : {type: 'histogram', interval: discretization, exclude_from_search, extended_bounds}
  }, [filter, fromDisplayUnit, isTime, discretization, nBins, autorange])
  const agg = useAgg(xAxis.quantity, visible && !disableHistogram, `${aggId}_histogram`, aggHistogramConfig)
  useEffect(() => {
    if (!isNil(agg)) {
      firstLoad.current = false
    }
  }, [agg])

  // Aggregation when the statistics are disabled: a simple min_max aggregation
  // is enough in order to get the slider range.
  const aggSliderConfig = useMemo(() => ({type: 'min_max', exclude_from_search: true}), [])
  const aggSlider = useAgg(xAxis.quantity, visible && disableHistogram, `${aggId}_slider`, aggSliderConfig)

  // Determine the global minimum and maximum values
  const [minGlobal, maxGlobal] = useMemo(() => {
    let minGlobal
    let maxGlobal
    if (disableHistogram) {
      minGlobal = aggSlider?.data?.[0]
      maxGlobal = aggSlider?.data?.[1]
    } else {
      const nBuckets = agg?.data?.length || 0
      if (nBuckets === 1) {
        minGlobal = agg.data[0].value
        maxGlobal = minGlobal
      } else if (nBuckets > 1) {
        for (const bucket of agg.data) {
          if (isNil(minGlobal)) {
            minGlobal = bucket.value
          }
          maxGlobal = bucket.value + (discretization ? 0 : agg.interval)
        }
        if (isNil(minGlobal)) {
          minGlobal = agg.data[0].value
        }
        if (isNil(maxGlobal)) {
          maxGlobal = agg.data[agg.data.length - 1].value + (discretization ? 0 : agg.interval)
        }
      }
    }
    firstRender.current = false
    return [minGlobal, maxGlobal]
  }, [agg, aggSlider, disableHistogram, discretization])

  const stepHistogram = agg?.interval
  const unavailable = isNil(minGlobal) || isNil(maxGlobal) || isNil(range)
  const disabled = unavailable

  // Determine the step value for sliders. Notice that this does not have to
  // match with the histogram binning, and that we want to do call the
  // getInterval function on the display unit range.
  const stepSlider = useMemo(() => {
    if (discretization) {
      return discretization
    }
    if (isNil(maxLocal) || isNil(minLocal)) {
      return undefined
    }
    const rangeSI = maxLocal - minLocal
    const range = new Quantity(rangeSI, unitStorage).to(xAxis.unit).value()
    const intervalCustom = getInterval(range, nSteps, xAxis.dtype)
    return new Quantity(intervalCustom, xAxis.unit).to(unitStorage).value()
  }, [maxLocal, minLocal, discretization, nSteps, xAxis.dtype, xAxis.unit, unitStorage])

  // When filter changes, the plot data should not be updated.
  useEffect(() => {
    if (autorange) loading.current = true
  }, [filter, autorange])

  // Once the aggregation data arrives, the plot data can be updated.
  useEffect(() => {
    if (!isNil(agg)) {
      loading.current = false
    }
  }, [agg])

  // We refresh the plot information only when all of the data is updated.
  // Essentially this syncs the axis limits and data so that they are updated
  // simultaneously.
  useEffect(() => {
    if (loading.current || isNil(agg?.data)) {
      return
    }

    setPlotData({
      xAxis: {
        quantity: xAxis.quantity,
        unit: xAxis.unit,
        unitStorage: unitStorage,
        dtype: xAxis.dtype,
        title: xAxis.title,
        min: minLocal,
        max: maxLocal
      },
      step: stepHistogram,
      data: agg.data
    })
  }, [loading, nBins, agg, minLocal, maxLocal, stepHistogram, unitStorage, xAxis.quantity, xAxis.unit, xAxis.dtype, xAxis.title])

  // Function for converting search values into the currently selected unit
  // system.
  const toInternal = useCallback(filter => {
    return (!isTime && unitStorage)
      ? formatNumber(new Quantity(filter, unitStorage).to(xAxis.unit).value())
      : filter
  }, [unitStorage, isTime, xAxis.unit])

  // If no filter has been specified by the user, the range is automatically
  // adjusted according to global min/max of the field. If filter is set, the
  // slider value is set according to it.
  useEffect(() => {
    let lte
    let gte
    let min
    let max
    let minInc = true
    let maxInc = true

    // Helper functions for determining the min/max boundaries.
    const limit = (global, filter, min) => {
      const func = min ? Math.min : Math.max
      if (autorange) {
        return filter
      } else {
        return func(global, filter)
      }
    }
    const limitMin = (global, filter) => limit(global, filter, true)
    const limitMax = (global, filter) => limit(global, filter, false)

    if (!isNil(minGlobal) && !isNil(maxGlobal)) {
      // When no filter is set, use the whole available range
      if (isNil(filter)) {
        gte = minGlobal
        lte = maxGlobal
        min = minGlobal
        max = maxGlobal
      // A single specific value is given
      } else if (filter instanceof Quantity) {
        gte = filter.to(unitStorage).value()
        lte = filter.to(unitStorage).value()
        min = limitMin(minGlobal, gte)
        max = limitMax(maxGlobal, lte)
      // A range is given. For visualization purposes open-ended queries are
      // displayed as well, although making such queries is currently not
      // supported.
      } else {
        minInc = isNil(filter.gt)
        maxInc = isNil(filter.lt)
        gte = filter.gte || filter.gt
        lte = filter.lte || filter.lt
        gte = gte instanceof Quantity ? gte.to(unitStorage).value() : gte
        lte = lte instanceof Quantity ? lte.to(unitStorage).value() : lte

        if (isNil(gte)) {
          min = limitMin(minGlobal, lte)
          gte = min
        } else {
          min = limitMin(minGlobal, gte)
        }
        if (isNil(lte)) {
          max = limitMax(maxGlobal, gte)
          lte = max
        } else {
          max = limitMax(maxGlobal, lte)
        }
      }
      minRef.current = min
      maxRef.current = max
      setMinInclusive(minInc)
      setMaxInclusive(maxInc)
      setMinLocal(min)
      setMaxLocal(max)
      setRange({gte: gte, lte: lte})
      setMinInput(toInternal(gte))
      setMaxInput(toInternal(lte))
      setMaxError()
      setMinError()
    }
  }, [minGlobal, maxGlobal, filter, unitStorage, toInternal, def, isTime, autorange])

  // Returns whether the given range is an acceptable value to be queried and
  // displayed.
  const validateRange = useCallback((range) => {
    if (isTime) {
      return range[1] - range[0] >= 60000
    }
    return true
  }, [isTime])

  // Handles changes in the min input
  const handleMinChange = useCallback((value) => {
    if (autorange) loading.current = true
    let val
    if (isTime) {
      value.setSeconds(0, 0)
      val = getTime(value)
    } else {
      val = value.target?.value
    }
    minInputChanged.current = true
    setMinInput(val)
    setMaxError()
    setMinError()
  }, [isTime, autorange])

  // Handles changes in the max input
  const handleMaxChange = useCallback((value) => {
    if (autorange) loading.current = true
    loading.current = true
    let val
    if (isTime) {
      value.setSeconds(0, 0)
      val = getTime(value)
    } else {
      val = value.target?.value
    }
    maxInputChanged.current = true
    setMaxInput(val)
    setMaxError()
    setMinError()
  }, [isTime, autorange])

  // Called when min values are submitted through the input field
  const handleMinSubmit = useCallback((value) => {
    let error
    if (isTime) {
      if (!isValid(value)) {
        error = 'Invalid start date.'
      } else {
        value = getTime(value)
      }
    } else {
      const number = Number.parseFloat(value)
      if (isNaN(number)) {
        error = 'Invalid minimum value.'
      } else {
        value = toStorageUnit(number)
      }
    }
    if ((isTime ? value : value.value) > range.lte) {
      error = 'Cannot be larger than maximum value.'
    }
    if (error) {
      setMinError(error)
    } else {
      minInputChanged.current = false
      setFilter(old => {
        return {
          gte: value,
          lte: fromStorageUnit(isNil(old?.lte) ? maxLocal : old.lte)
        }
      })
    }
  }, [isTime, setFilter, fromStorageUnit, toStorageUnit, maxLocal, range])

  // Called when max values are submitted through the input field
  const handleMaxSubmit = useCallback((value) => {
    let error
    if (isTime) {
      if (!isValid(value)) {
        error = 'Invalid end date.'
      } else {
        value = getTime(value)
      }
    } else {
      const number = Number.parseFloat(value)
      if (isNaN(number)) {
        error = 'Invalid maximum value.'
      } else {
        value = toStorageUnit(number)
      }
    }
    if ((isTime ? value : value.value) < range.gte) {
      error = 'Cannot be smaller than minimum value.'
    }
    if (error) {
      setMaxError(error)
    } else {
      maxInputChanged.current = false
      setFilter(old => {
        return {
          gte: fromStorageUnit(isNil(old?.gte) ? minLocal : old.gte),
          lte: value
        }
      })
    }
  }, [isTime, setFilter, fromStorageUnit, toStorageUnit, minLocal, range])

  // Handle range commit: Set the filter when mouse is released on a slider.
  // Notice that we cannot rely on the value given by the slider event: it may
  // correspond to a non-valid range. Instead we need to use the latest valid
  // value saved when the value has changed.
  const handleRangeCommit = useCallback((event) => {
    const value = validRange.current
    if (!isNil(value) && rangeChanged.current) {
      setFilter({
        gte: fromStorageUnit(value[0]),
        lte: fromStorageUnit(value[1])
      })
      rangeChanged.current = false
    }
  }, [setFilter, fromStorageUnit])

  // Handle range change: only change the rendered values, filter is send with
  // handleRangeCommit.
  const handleRangeChange = useCallback((event, value, validate = true) => {
    if (autorange) loading.current = true
    const valid = !validate || validateRange(value)
    if (valid) {
      rangeChanged.current = true
      validRange.current = value
      setRange({gte: value[0], lte: value[1]})
      setMinInput(toInternal(value[0]))
      setMaxInput(toInternal(value[1]))
      setMinInclusive(true)
      setMaxInclusive(true)
      setMaxError()
      setMinError()
    }
  }, [validateRange, validRange, toInternal, autorange])

  // Handle min input field blur: only if the input has changed, the filter will
  // be submitted.
  const handleMinBlur = useCallback(() => {
    minInputChanged.current && handleMinSubmit(minInput)
  }, [minInput, handleMinSubmit, minInputChanged])

  // Handle max input field blur: only if the input has changed, the filter will
  // be submitted.
  const handleMaxBlur = useCallback(() => {
    maxInputChanged.current && handleMaxSubmit(maxInput)
  }, [maxInput, handleMaxSubmit, maxInputChanged])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <PlotHistogram
      bins={plotData?.data}
      xAxis={plotData?.xAxis}
      step={plotData?.step}
      minXInclusive={minInclusive}
      maxXInclusive={maxInclusive}
      disabled={disabled}
      scale={scale}
      nBins={nBins}
      range={range}
      highlight={highlight}
      discretization={discretization}
      tooltipLabel={repeats ? 'value' : 'entry'}
      dtypeY={DType.Int}
      onRangeChange={handleRangeChange}
      onRangeCommit={handleRangeCommit}
      onMinBlur={handleMinBlur}
      onMaxBlur={handleMaxBlur}
      onMinChange={handleMinChange}
      onMaxChange={handleMaxChange}
      onMinSubmit={handleMinSubmit}
      onMaxSubmit={handleMaxSubmit}
      classes={{histogram: classes?.histogram}}
      onClick={(event, value) => {
        handleRangeChange(event, value, false)
        handleRangeCommit(event, value)
      }}
      minError={minError}
      maxError={maxError}
      minInput={minInput}
      maxInput={maxInput}
      minLocal={minLocal}
      maxLocal={maxLocal}
      showinput={showinput}
      stepSlider={stepSlider}
      disableHistogram={disableHistogram}
      disableXTitle={disableXTitle}
      data-testid={`${testID}-histogram`}
    />
  </div>
})

Range.propTypes = {
  xAxis: PropTypes.object,
  /* Target number of steps for the slider that is shown when statistics are
   * disabled. The actual number may vary, as the step is chosen to be a
   * human-readable value that depends on the range and the unit. */
  nSteps: PropTypes.number,
  /* Number of bins for the histogram (and the slider) when statistics are
   * enabled. */
  nBins: PropTypes.number,
  visible: PropTypes.bool,
  /* The statistics scaling */
  scale: PropTypes.string,
  /* Whether the histogram is disabled */
  disableHistogram: PropTypes.bool,
  /* Whether the x title is disabled */
  disableXTitle: PropTypes.bool,
  /* Set the range automatically according to data. */
  autorange: PropTypes.bool,
  /* Show the input fields for min and max value */
  showinput: PropTypes.bool,
  aggId: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

Range.defaultProps = {
  nSteps: 20,
  nBins: 30,
  aggId: 'default',
  'data-testid': 'inputrange'
}

/**
 * A small wrapper around Range for use in the filter menus.
 */
const useInputRangeStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'column'
  },
  histogram: {
    height: '8rem'
  }
}))
const InputRange = React.memo(({
  label,
  quantity,
  description,
  nSteps,
  visible,
  initialScale,
  nBins,
  disableHistogram,
  initialAutorange,
  aggId,
  className,
  'data-testid': testID
}) => {
  const {filterData} = useSearchContext()
  const {units} = useUnitContext()
  const styles = useInputRangeStyles()
  const [scale, setScale] = useState(initialScale || filterData[quantity].scale)
  const dtype = filterData[quantity].dtype
  const isTime = dtype === DType.Timestamp
  const [autorange, setAutorange] = useState(isNil(initialAutorange) ? isTime : initialAutorange)
  const x = useMemo(() => (
    {
      quantity,
      dtype,
      unit: new Unit(filterData[quantity]?.unit || 'dimensionless').toSystem(units)
    }
  ), [quantity, filterData, dtype, units])

  // Determine the description and title
  const def = filterData[quantity]
  const descFinal = description || def?.description || ''
  const labelFinal = label || def?.label

  // Component for enabling/disabling autorange
  const actions = <ActionCheckbox
    value={autorange}
    label="zoom"
    tooltip={autorangeDescription}
    onChange={(value) => setAutorange(value)}
  />

  return <div className={clsx(styles.root, className)}>
    <InputHeader
      label={labelFinal}
      quantity={quantity}
      description={descFinal}
      scale={scale}
      onChangeScale={setScale}
      disableStatistics={disableHistogram}
      actions={actions}
    />
    <Range
      xAxis={x}
      nSteps={nSteps}
      visible={visible}
      scale={scale}
      nBins={nBins}
      disableHistogram={disableHistogram}
      disableXTitle
      autorange={autorange}
      showinput
      aggId={aggId}
      classes={{histogram: styles.histogram}}
      data-testid={testID}
    />
  </div>
})

InputRange.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  nSteps: PropTypes.number,
  nBins: PropTypes.number,
  visible: PropTypes.bool,
  initialScale: PropTypes.string,
  disableHistogram: PropTypes.bool,
  initialAutorange: PropTypes.bool,
  aggId: PropTypes.string,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

InputRange.defaultProps = {
  nSteps: 20,
  nBins: 30,
  aggId: 'default',
  'data-testid': 'inputrange'
}

export default InputRange
