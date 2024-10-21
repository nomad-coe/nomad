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
import React, { useRef, useMemo, useCallback } from 'react'
import clsx from 'clsx'
import { range as rangeLodash, isNil, clamp } from 'lodash'
import { useRecoilValue } from 'recoil'
import { Slider, InputAdornment } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import { KeyboardDateTimePicker } from '@material-ui/pickers'
import { pluralize, formatInteger, DType } from '../../utils'
import { Quantity } from '../units/Quantity'
import InputUnavailable from '../search/input/InputUnavailable'
import { InputTextField } from '../search/input/InputText'
import Placeholder from '../visualization/Placeholder'
import PlotAxis from './PlotAxis'
import PlotBar from './PlotBar'
import FilterTitle from '../search/FilterTitle'
import { guiState } from '../GUIMenu'
import PropTypes from 'prop-types'
import { getScaler } from './common'
import { dateFormat } from '../../config'

/**
 * An interactive histogram for numeric values.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    boxSizing: 'border-box'
  },
  container: {
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center'
  },
  histogram: {
    width: '100%',
    height: '100%'
  },
  overflow: {
    flex: '1 1 auto',
    width: '100%',
    height: '100%'
  },
  grid: {
    width: '100%',
    height: '100%',
    boxSizing: 'border-box',
    display: 'grid',
    gridTemplateColumns: 'auto 1fr',
    gridTemplateRows: '1fr auto',
    paddingRight: theme.spacing(0.75),
    paddingBottom: theme.spacing(0.25)
  },
  plot: {
    gridColumn: 2,
    gridRow: 1,
    position: 'relative'
  },
  canvas: {
    width: '100%',
    height: '100%',
    position: 'relative',
    overflow: 'hidden'
  },
  slider: {
    zIndex: 10,
    position: 'absolute',
    bottom: -14,
    left: 0,
    right: 0
  },
  yaxis: {
    gridColumn: 1,
    gridRow: 1,
    height: '100%'
  },
  square: {
    gridColumn: 1,
    gridRow: 2
  },
  xaxis: {
    gridColumn: 2,
    gridRow: 2
  },
  inputField: {
    marginTop: 0,
    marginBottom: 0,
    flexGrow: 1,
    minWidth: '6.1rem',
    maxWidth: '8.5rem'
  },
  inputFieldDate: {
    maxWidth: '15rem'
  },
  thumb: {
    '&:active': {
      boxShadow: '0px 0px 0px 12px rgb(0, 141, 195, 16%)'
    },
    '&:focusVisible': {
      boxShadow: '0px 0px 0px 6px rgb(0, 141, 195, 16%)'
    }
  },
  title: {
    flexGrow: 1,
    marginLeft: theme.spacing(0.5),
    marginRight: theme.spacing(0.5)
  },
  // The following two styles are needed in order for the TextInput to
  // transition from having a label to not having a label. The MUI component
  // does not otherwise transition correctly.
  inputMargin: {
    paddingTop: 10,
    paddingBottom: 11
  },
  adornment: {
    marginTop: '0px !important'
  },
  spacer: {
    height: '3rem',
    flex: '1 1 100%',
    paddingLeft: '18px',
    paddingRight: '18px',
    display: 'flex',
    alignItems: 'center'
  },
  row: {
    marginTop: theme.spacing(0.25),
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    justifyContent: 'center',
    alignItems: 'flex-start'
  },
  axisTitle: {
    fontSize: '0.75rem'
  }
}))
const PlotHistogram = React.memo(({
  xAxis,
  yAxis,
  bins,
  range,
  step,
  nBins,
  discretization,
  dtypeY,
  disabled,
  tooltipLabel,
  highlight,
  disableSlider,
  onRangeChange,
  onRangeCommit,
  onMinChange,
  onMaxChange,
  onMinBlur,
  onMaxBlur,
  onMinSubmit,
  onMaxSubmit,
  onClick,
  minXInclusive,
  maxXInclusive,
  className,
  classes,
  showInput,
  minError,
  maxError,
  minInput,
  maxInput,
  minLocal,
  maxLocal,
  stepSlider,
  disableHistogram,
  disableXTitle,
  'data-testid': testID
}) => {
  const styles = useStyles(classes)
  const titleClasses = {text: styles.axisTitle}
  const useDynamicStyles = makeStyles((theme) => {
    const color = highlight ? theme.palette.secondary.main : theme.palette.primary.main
    return {
      thumb: {
        'border': `2px solid ${color}`,
        '&[data-index="0"]': {
          backgroundColor: minXInclusive ? color : 'white'
        },
        '&[data-index="1"]': {
          backgroundColor: maxXInclusive ? color : 'white'
        }
      }
    }
  })
  const dynamicStyles = useDynamicStyles()

  const aggIndicator = useRecoilValue(guiState('aggIndicator'))
  const oldRangeRef = useRef()
  const artificialRange = 1
  const isArtificial = !step && bins?.length === 1
  const isTime = xAxis.dtype === DType.Timestamp
  step = isArtificial ? artificialRange / nBins : step

  // Determine the final maximum of x-axis. If the values are discrete, the
  // maximum x value needs to be increased.
  let maxX = isArtificial ? xAxis.min + artificialRange : xAxis.max
  if (discretization) {
    if (!isNil(maxX)) {
      maxX += step
    }
  }

  // Determine the final bin data and the y axis limits.
  const [finalBins, minY, maxY] = useMemo(() => {
    if (isNil(bins)) {
      return [null, null, null]
    }
    const minY = 0
    const maxY = Math.max(...bins.map(item => item.count))
    const scaler = getScaler(yAxis.scale, [minY, maxY])
    const finalBins = bins.map((bucket) => {
      return {
        ...bucket,
        start: bucket.value,
        end: bucket.value + step,
        scale: scaler(bucket.count)
      }
    })
    return [finalBins, minY, maxY]
  }, [bins, step, yAxis.scale])

  // Transforms the original range into an internal range used for
  // visualization.
  const toInternalRange = useCallback((range) => {
    return (discretization && !isNil(range?.[1]))
      ? [range[0], range[1] + step]
      : range
  }, [discretization, step])

  // Transforms the internal range into the original range.
  const toExternalRange = useCallback((range) => {
    return (discretization && !isNil(range?.[1]))
      ? [range[0], range[1] - step]
      : range
  }, [discretization, step])

  const rangeInternal = useMemo(() => {
    return toInternalRange([range?.gte, range?.lte])
  }, [range, toInternalRange])

  // Turns the internal range values back to the external counterparts before
  // sending the event.
  const handleRangeCommit = useCallback((event, range) => {
    onRangeCommit && onRangeCommit(event, toExternalRange(range))
  }, [onRangeCommit, toExternalRange])

  // Handle clicks on the histogram bins
  const handleClick = useCallback((event, item) => {
    onClick && onClick(event, isArtificial
      ? [item.start, item.start]
      : [item.start, discretization ? item.start : item.end])
  }, [isArtificial, onClick, discretization])

  // Turns the internal range values back to the external counterparts before
  // sending the event.
  const handleRangeChange = useCallback((event, range) => {
    if (!onRangeChange) {
      return
    }
    range = toExternalRange(range)
    // Check if the value has actually changed: MUI Slider does not do this
    // automatically.
    if (oldRangeRef.current &&
      oldRangeRef.current[0] === range[0] &&
      oldRangeRef.current[1] === range[1]) {
      return
    }
    oldRangeRef.current = range
    onRangeChange(event, range)
  }, [onRangeChange, toExternalRange])

  // Calculates the selection value for a bin.
  const calculateSelection = useCallback((bin, range) => {
    if (!highlight) {
      return false
    }
    if (isArtificial) {
      return bin.start === rangeInternal[0]
    }
    const interval = bin.end - bin.start
    const min = clamp((range[0] - bin.start) / interval, 0, 1)
    const max = clamp((range[1] - bin.start) / interval, 0, 1)
    const size = max - min
    if (size <= 0) {
      return false
    }
    if (min <= 0 && max >= 1) {
      return true
    }
    return [min, max]
  }, [highlight, isArtificial, rangeInternal])

  // Create the plot once items are ready
  const plot = useMemo(() => {
    if (isNil(finalBins) || isNil(xAxis.min) || isNil(maxX)) {
      return null
    }
    return <div className={styles.canvas}>
      {Object.values(finalBins).map((item, i) => {
        return <PlotBar
          startX={(item.start - xAxis.min) / (maxX - xAxis.min)}
          endX={(item.end - xAxis.min) / (maxX - xAxis.min)}
          startY={0}
          endY={item.scale}
          key={i}
          tooltip={`${formatInteger(item.count)} ${pluralize(tooltipLabel, item.count, false)}`}
          vertical
          disableValue
          selected={calculateSelection(item, rangeInternal)}
          onClick={(event) => handleClick(event, item)}
        />
      })}
    </div>
  }, [finalBins, xAxis.min, maxX, styles.canvas, calculateSelection, rangeInternal, handleClick, tooltipLabel])

  // Create x-axis once axis range is ready
  const xaxis = useMemo(() => {
    if (isNil(finalBins) || isNil(xAxis.min) || isNil(maxX)) {
      return null
    }
    let labels

    // One bin is shown for the artificial values
    if (isArtificial) {
      const offset = discretization ? 0.5 * step : 0
      labels = [{
        label: finalBins[0].start,
        pos: finalBins[0].start + offset
      }]
    // Discrete values get label at the center of the bin.
    } else if (discretization) {
      const start = step * Math.ceil(xAxis.min / step)
      const end = step * Math.floor(maxX / step)
      labels = rangeLodash(start, end).map(x => ({
        label: x,
        pos: x + 0.5 * step
      }))
    }

    const min = new Quantity(xAxis.min, xAxis.unitStorage).to(xAxis.unit).value()
    const max = new Quantity(maxX, xAxis.unitStorage).to(xAxis.unit).value()

    return <PlotAxis
      placement="bottom"
      min={min}
      max={max}
      mode="scientific"
      labels={labels}
      labelWidth={45}
      overflowLeft={25}
      dtype={xAxis.dtype}
      className={styles.xaxis}
    />
  }, [finalBins, xAxis.min, xAxis.unitStorage, maxX, discretization, isArtificial, xAxis.dtype, styles.xaxis, xAxis.unit, step])

  const [minAdornment, maxAdornment] = useMemo(() => {
    return disableHistogram
      ? [undefined, undefined]
      : [
      {
        startAdornment: <InputAdornment
          position="start"
          classes={{positionStart: styles.adornment}}
        >min:</InputAdornment>,
        classes: {inputMarginDense: styles.inputMargin}
      },
      {
        startAdornment: <InputAdornment
          position="start"
          classes={{positionStart: styles.adornment}}
        >max:</InputAdornment>,
        classes: {inputMarginDense: styles.inputMargin}
      }
    ]
  }, [disableHistogram, styles])

  // Determine the min input component
  let inputMinField
  if (xAxis.dtype === DType.Timestamp) {
    inputMinField = <KeyboardDateTimePicker
      error={!!minError}
      disabled={disabled}
      helperText={minError}
      ampm={false}
      className={clsx(styles.inputField, styles.inputFieldDate)}
      variant="inline"
      inputVariant="filled"
      label="Start time"
      format={`${dateFormat} kk:mm`}
      value={minInput}
      invalidDateMessage=""
      InputAdornmentProps={{ position: 'end' }}
      onAccept={(date) => {
        onMinChange(date)
        onMinSubmit(date)
      }}
      onChange={onMinChange}
      onBlur={onMinBlur}
      onKeyDown={(event) => { if (event.key === 'Enter') { onMinSubmit(minInput) } }}
    />
  } else {
    inputMinField = <InputTextField
      error={!!minError}
      disabled={disabled}
      helperText={minError}
      className={styles.inputField}
      InputProps={minAdornment}
      label={disableHistogram ? 'min' : ' '}
      value={minInput}
      onChange={onMinChange}
      onBlur={onMinBlur}
      onKeyDown={(event) => { if (event.key === 'Enter') { onMinSubmit(minInput) } }}
    />
  }

  // Determine the max input component
  let inputMaxField
  if (xAxis.dtype === DType.Timestamp) {
    inputMaxField = <KeyboardDateTimePicker
      error={!!maxError}
      disabled={disabled}
      helperText={maxError}
      ampm={false}
      className={clsx(styles.inputField, styles.inputFieldDate)}
      variant="inline"
      inputVariant="filled"
      label="End time"
      format={`${dateFormat} kk:mm`}
      value={maxInput}
      invalidDateMessage=""
      InputAdornmentProps={{ position: 'end' }}
      onAccept={(date) => {
        onMaxChange(date)
        onMaxSubmit(date)
      }}
      onChange={onMaxChange}
      onBlur={onMaxBlur}
      onKeyDown={(event) => { if (event.key === 'Enter') { onMaxSubmit(maxInput) } }}
    />
  } else {
    inputMaxField = <InputTextField
      error={!!maxError}
      disabled={disabled}
      helperText={maxError}
      className={styles.inputField}
      InputProps={maxAdornment}
      label={disableHistogram ? 'max' : ' '}
      value={maxInput}
      onChange={onMaxChange}
      onBlur={onMaxBlur}
      onKeyDown={(event) => { if (event.key === 'Enter') { onMaxSubmit(maxInput) } }}
    />
  }

  // Create y-axis once axis range is ready.
  const yaxis = useMemo(() => {
    if (isNil(minY) || isNil(maxY)) {
      return null
    }
    return <PlotAxis
      placement="left"
      min={minY}
      max={maxY}
      mode='SI'
      scale={yAxis.scale}
      dtype={dtypeY}
      className={styles.yaxis}
    />
  }, [dtypeY, maxY, minY, yAxis.scale, styles.yaxis])

  // Determine the final component to show.
  let histComp
  if (bins?.length === 0) {
    histComp = <InputUnavailable/>
  } else if (isNil(bins) && aggIndicator === 'on') {
    histComp = <Placeholder
      variant="rect"
      data-testid={`${testID}-placeholder`}
      margin={0}
    />
  } else {
    histComp = <div className={styles.overflow}>
      <div className={styles.grid}>
        {yaxis}
        <div className={styles.plot}>
          {plot}
          {!disableSlider &&
            <Slider
              disabled={disabled || isArtificial || (maxX - xAxis.min) === discretization}
              color={highlight ? 'secondary' : 'primary'}
              min={xAxis.min}
              max={maxX}
              step={step}
              value={rangeInternal}
              onChange={handleRangeChange}
              onChangeCommitted={handleRangeCommit}
              valueLabelDisplay="off"
              classes={{thumb: styles.thumb && dynamicStyles.thumb}}
              className={styles.slider}
            />
          }
        </div>
        <div className={styles.square} />
        {xaxis}
      </div>
  </div>
  }

  const titleComp = <div className={styles.title}>
    <FilterTitle
      variant="subtitle2"
      classes={titleClasses}
      quantity={xAxis.quantity}
      label={xAxis.title}
      unit={xAxis.unit}
      noWrap={false}
    />
  </div>

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div className={styles.container}>
      {!disableHistogram && <div className={clsx(styles.histogram, classes?.histogram)}>{histComp}</div>}
      {!disableXTitle && titleComp}
      <div className={styles.row}>
        {showInput
          ? <>
            {inputMinField}
            {(disableHistogram && !isTime)
              ? <div className={styles.spacer}>
                  <Slider
                    disabled={disabled || (minLocal === maxLocal)}
                    color={highlight ? 'secondary' : 'primary'}
                    min={minLocal}
                    max={maxLocal}
                    step={stepSlider}
                    value={[range.gte, range.lte]}
                    onChange={onRangeChange}
                    onChangeCommitted={onRangeCommit}
                    valueLabelDisplay="off"
                  />
                </div>
              : <div className={styles.title} />
            }
            {inputMaxField}
            </>
          : null
        }
      </div>
    </div>
  </div>
})

PlotHistogram.propTypes = {
  xAxis: PropTypes.object,
  yAxis: PropTypes.object,
  /* The bins data to show. */
  bins: PropTypes.arrayOf(PropTypes.shape({
    value: PropTypes.number,
    count: PropTypes.number
  })),
  range: PropTypes.object,
  step: PropTypes.number,
  /* The number of bins: required only when showing a dummy bin for histograms
   * that have no width. */
  nBins: PropTypes.number,
  /* Discretization of the values. */
  discretization: PropTypes.number,
  dtypeY: PropTypes.string,
  disabled: PropTypes.bool,
  /* The label to show for the tooltips */
  tooltipLabel: PropTypes.string,
  /* Whether the selected region should be highlighted */
  highlight: PropTypes.bool,
  /* Disable the intervative slider shown on top of the plot */
  disableSlider: PropTypes.bool,
  /* Function to call when a histogram bar has been clicked */
  onClick: PropTypes.func,
  /* Whether the min slider is inclusive (=min value is included) */
  minXInclusive: PropTypes.bool,
  /* Whether the max slider is inclusive (=max value is included) */
  maxXInclusive: PropTypes.bool,
  showInput: PropTypes.bool,
  onRangeChange: PropTypes.func,
  onRangeCommit: PropTypes.func,
  onMinChange: PropTypes.func,
  onMaxChange: PropTypes.func,
  onMinSubmit: PropTypes.func,
  onMaxSubmit: PropTypes.func,
  maxError: PropTypes.any,
  minError: PropTypes.any,
  maxInput: PropTypes.any,
  minInput: PropTypes.any,
  maxLocal: PropTypes.any,
  minLocal: PropTypes.any,
  stepSlider: PropTypes.any,
  disableHistogram: PropTypes.bool,
  disableXTitle: PropTypes.bool,
  onMinBlur: PropTypes.func,
  onMaxBlur: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default PlotHistogram
