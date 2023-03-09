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
import { Slider } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import { pluralize, formatInteger } from '../../utils'
import { Unit } from '../../units'
import InputUnavailable from '../search/input/InputUnavailable'
import Placeholder from '../visualization/Placeholder'
import PlotAxis from './PlotAxis'
import PlotBar from './PlotBar'
import { guiState } from '../GUIMenu'
import PropTypes from 'prop-types'
import { getScaler } from './common'

/**
 * An interactive histogram for numeric values.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    boxSizing: 'border-box'
  },
  overflow: {
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
    paddingBottom: theme.spacing(0.5)
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
  thumb: {
    '&:active': {
      boxShadow: '0px 0px 0px 12px rgb(0, 141, 195, 16%)'
    },
    '&:focusVisible': {
      boxShadow: '0px 0px 0px 6px rgb(0, 141, 195, 16%)'
    }
  }
}))
const PlotHistogram = React.memo(({
  bins,
  range,
  minX,
  maxX,
  unitX,
  step,
  nBins,
  scale,
  discretization,
  dtypeX,
  dtypeY,
  disabled,
  tooltipLabel,
  highlight,
  disableSlider,
  onRangeChange,
  onRangeCommit,
  onClick,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = useStyles(classes)
  const scaler = useMemo(() => getScaler(scale), [scale])
  const aggIndicator = useRecoilValue(guiState('aggIndicator'))
  const oldRangeRef = useRef()
  const artificialRange = 1
  const isArtificial = !step && bins?.length === 1
  step = isArtificial ? artificialRange / nBins : step
  maxX = isArtificial ? minX + artificialRange : maxX

  // Determine the final bin data and the y axis limits.
  const [finalBins, minY, maxY] = useMemo(() => {
    if (isNil(bins)) {
      return [null, null, null]
    }
    const minY = 0
    const maxY = Math.max(...bins.map(item => item.count))
    const finalBins = bins.map((bucket) => {
      return {
        ...bucket,
        start: bucket.value,
        end: bucket.value + step,
        scale: scaler(bucket.count / maxY) || 0
      }
    })
    return [finalBins, minY, maxY]
  }, [bins, scaler, step])

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

  // If the values are discrete, the maximum x value needs to be increased.
  const rangeInternal = useMemo(() => {
    return toInternalRange([range?.gte, range?.lte])
  }, [range, toInternalRange])
  if (discretization) {
    if (!isNil(maxX)) {
      maxX += step
    }
  }

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
    if (isNil(finalBins) || isNil(minX) || isNil(maxX)) {
      return null
    }
    return <div className={styles.canvas}>
      {Object.values(finalBins).map((item, i) => {
        return <PlotBar
          startX={(item.start - minX) / (maxX - minX)}
          endX={(item.end - minX) / (maxX - minX)}
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
  }, [finalBins, minX, maxX, styles.canvas, calculateSelection, rangeInternal, handleClick, tooltipLabel])

  // Create x-axis once axis range is ready
  const xaxis = useMemo(() => {
    if (isNil(finalBins) || isNil(minX) || isNil(maxX)) {
      return null
    }

    let labels
    // Automatic labelling is used for continuous values
    if (!discretization && !isArtificial) {
      labels = 10
    // One bin is shown for the artificial values
    } else if (isArtificial) {
      const offset = discretization ? 0.5 * step : 0
      labels = [{
        label: finalBins[0].start,
        pos: finalBins[0].start + offset
      }]
    // Discrete values get label at the center of the bin.
    } else {
      const start = step * Math.ceil(minX / step)
      const end = step * Math.floor(maxX / step)
      labels = rangeLodash(start, end).map(x => ({
        label: x,
        pos: x + 0.5 * step
      }))
    }

    return <PlotAxis
      placement="bottom"
      min={minX}
      max={maxX}
      unit={unitX}
      mode="scientific"
      labels={labels}
      labelWidth={45}
      overflowLeft={25}
      dtype={dtypeX}
      className={styles.xaxis}
    />
  }, [finalBins, minX, maxX, discretization, isArtificial, unitX, dtypeX, styles.xaxis, step])

  // Create y-axis once axis range is ready.
  const yaxis = useMemo(() => {
    if (isNil(minY) || isNil(maxY)) {
      return null
    }
    return <PlotAxis
      placement="left"
      min={minY}
      max={maxY}
      unit={new Unit('dimensionless')}
      mode='SI'
      labels={5}
      scale={scale}
      dtype={dtypeY}
      className={styles.yaxis}
    />
  }, [dtypeY, maxY, minY, scale, styles.yaxis])

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
              disabled={disabled || isArtificial || (maxX - minX) === discretization}
              color={highlight ? 'secondary' : 'primary'}
              min={minX}
              max={maxX}
              step={step}
              value={rangeInternal}
              onChange={handleRangeChange}
              onChangeCommitted={handleRangeCommit}
              valueLabelDisplay="off"
              classes={{thumb: styles.thumb}}
              className={styles.slider}
            />
          }
        </div>
        <div className={styles.square} />
        {xaxis}
      </div>
    </div>
  }

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    {histComp}
  </div>
})

PlotHistogram.propTypes = {
  /* The bins data to show. */
  bins: PropTypes.arrayOf(PropTypes.shape({
    value: PropTypes.number,
    count: PropTypes.number
  })),
  range: PropTypes.object,
  minX: PropTypes.number,
  maxX: PropTypes.number,
  unitX: PropTypes.any,
  step: PropTypes.number,
  /* The number of bins: required only when showing a dummy bin for histograms
   * that have no width. */
  nBins: PropTypes.number,
  /* Discretization of the values. */
  discretization: PropTypes.number,
  dtypeX: PropTypes.string,
  dtypeY: PropTypes.string,
  scale: PropTypes.string,
  disabled: PropTypes.bool,
  /* The label to show for the tooltips */
  tooltipLabel: PropTypes.string,
  /* Whether the selected region should be highlighted */
  highlight: PropTypes.bool,
  /* Disable the intervative slider shown on top of the plot */
  disableSlider: PropTypes.bool,
  /* Function to call when a histogram bar has been clicked */
  onClick: PropTypes.func,
  onRangeChange: PropTypes.func,
  onRangeCommit: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default PlotHistogram
