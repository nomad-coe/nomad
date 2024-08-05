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
import clsx from 'clsx'
import { format } from 'date-fns'
import { makeStyles } from '@material-ui/core/styles'
import { isArray, isNil } from 'lodash'
import { useResizeDetector } from 'react-resize-detector'
import { getScaler, getTicks } from './common'
import { formatNumber, DType } from '../../utils'
import PlotLabel from './PlotLabel'
import PlotTick from './PlotTick'
import PropTypes from 'prop-types'

/**
 * Represents a plot axis. Can be oriented to fit different sides of a plot
 * using the placement argument.
 */
const usePlotAxisStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    position: 'relative',
    overflow: 'hidden'
  },
  left: {
    flexDirection: 'row',
    height: '100%'
  },
  bottom: {
    flexDirection: 'column-reverse',
    width: '100%'
  },
  labelContainer: {
    position: 'relative'
  },
  labelPositioner: {
    position: 'absolute'
  },
  topTick: {
    position: 'absolute',
    top: 0,
    right: 0
  },
  verticalAxis: {
    position: 'absolute',
    width: 1,
    right: 0,
    backgroundColor: theme.palette.grey[500]
  },
  horizontalAxis: {
    position: 'absolute',
    height: 1,
    top: 0,
    backgroundColor: theme.palette.grey[500]
  }
}))

const PlotAxis = React.memo(({
  min,
  max,
  dtype,
  mode,
  decimals,
  scale,
  labels,
  labelHeight,
  labelWidth,
  placement,
  tickLength,
  labelPadding,
  overflowLeft,
  overflowRight,
  overflowTop,
  overflowBottom,
  className,
  classes,
  'data-testid': testID}) => {
  const styles = usePlotAxisStyles(classes)
  const {height, width, ref} = useResizeDetector()
  const orientation = {
    left: 'vertical',
    bottom: 'horizontal'
  }[placement]
  const axisSize = orientation === 'vertical'
    ? height
    : width
  const labelSize = (orientation === 'vertical'
    ? labelHeight
    : labelWidth)

  // Determine the correct scaler
  const scaler = useMemo(
    () => getScaler(scale, [min, max], [0, axisSize])
  , [scale, min, max, axisSize])

  // Determine styles that depend on overflow values
  overflowBottom = isNil(overflowBottom) ? 8 : overflowBottom
  overflowTop = isNil(overflowTop) ? (placement === 'bottom' ? 0 : 8) : overflowTop
  overflowLeft = isNil(overflowLeft) ? (placement === 'bottom' ? 0 : 8) : overflowLeft
  overflowRight = isNil(overflowRight) ? (placement === 'left' ? 0 : 8) : overflowRight
  const rootStyle = useMemo(() => {
    return {
      marginTop: -overflowTop,
      paddingTop: overflowTop,
      marginBottom: -overflowBottom,
      paddingBottom: overflowBottom,
      marginRight: -overflowRight,
      paddingRight: overflowRight,
      marginLeft: -overflowLeft,
      paddingLeft: overflowLeft
    }
  }, [overflowBottom, overflowLeft, overflowRight, overflowTop])
  const axisStyle = useMemo(() => {
    let axisStyle
    if (placement === 'left') {
      axisStyle = {
        top: overflowTop,
        bottom: overflowBottom - 1
      }
    } else if (placement === 'bottom') {
      axisStyle = {
        left: overflowLeft - 1,
        right: overflowRight - 1
      }
    }
    return axisStyle
  }, [overflowBottom, overflowLeft, overflowRight, overflowTop, placement])

  const labelList = useMemo(() => {
    if (!axisSize) {
      return []
    }

    // Function for formatting ticks
    const formatTick = (value) => {
      return dtype === DType.Timestamp
        ? format(value, 'MMM d')
        : formatNumber(value, dtype, mode, decimals)
    }

    // Manual ticks
    if (isArray(labels)) {
      return labels.map((tick) => ({...tick, label: formatTick(tick.label), pos: scaler(tick.pos) / axisSize}))
    }

    // On linear and log axes, the labels are spaced evenly, and the number of
    // labels is calculated from the available space. On non-linearly spaced
    // axes, we calculate the number of labels with formula:
    // axisSize - scaler(max - (max-min) / nLabels) = labelSize
    // -> nLabels = (max-min) / (max - scaler.invert(axisSize - labelSize)
    const padding = 10
    let nItems = (scale === 'linear' || scale === 'log')
      ? Math.floor(axisSize / (labelSize + padding))
      : Math.floor((max - min) / (max - scaler.invert(axisSize - (labelSize + padding))))

    // At least two labels should be attempted to be shown
    nItems = Math.max(2, nItems)

    // If the scale length is zero, show only one tick
    if (min === max) {
      return [{
        label: formatTick(min),
        pos: 0
      }]
    }

    // Get reasonable, human-readable ticks. the .ticks function from d3-scale
    // does not guarantee an upper limit to the number of ticks, so it cannot be
    // directly used.
    return getTicks(min, max, nItems, scale, dtype, mode, decimals)
      .map(({tick, value}) => {
        return {
          label: tick,
          pos: scaler(value) / axisSize
        }
      })
  }, [axisSize, dtype, labelSize, labels, max, min, scaler, mode, decimals, scale])

  // Here we estimate the maximum label width. This is a relatively simple
  // approximattion calculated using the font size. A more reliable way would to
  // to actually render to label and measure it's size on the DOM, but this is
  // more expensive.
  const finalLabelWidth = useMemo(() => {
    if (labelWidth) {
      return labelWidth
    }
    const maxCharacters = labelList?.length
      ? Math.max(...labelList.map(label => label.label.length))
      : 0
    const width = tickLength + labelPadding + maxCharacters * 0.56 * labelHeight
    return width
  }, [labelWidth, labelList, tickLength, labelPadding, labelHeight])

  return <div
    ref={ref}
    className={clsx(className, styles[placement], styles.root)}
    data-testid={testID || `plot-axis-${placement}`}
    style={rootStyle}
  >
    <div
      className={styles.labelContainer}
      style={orientation === 'vertical'
        ? {width: finalLabelWidth, height: '100%'}
        : {height: labelHeight + labelPadding + tickLength, width: '100%'}
      }
    >
      {labelList.map(({pos, label}, i) => <div
        key={i}
        className={styles.labelPositioner}
        style={orientation === 'vertical'
          ? {bottom: `${pos * 100}%`, right: 0}
          : {left: `${pos * 100}%`, top: 0}
        }
      >
        <PlotLabel
          label={label}
          labelPadding={labelPadding}
          size={labelHeight}
          tickLength={tickLength}
          orientation={orientation}
        />
      </div>)}
      {orientation === 'vertical' && <div className={styles.topTick}>
        <PlotTick orientation={orientation} length={tickLength}/>
      </div>}
    </div>
    <div className={styles[`${orientation}Axis`]} style={axisStyle}/>
  </div>
})

PlotAxis.propTypes = {
  min: PropTypes.number,
  max: PropTypes.number,
  scale: PropTypes.string,
  mode: PropTypes.oneOf(['scientific', 'SI', 'standard']),
  decimals: PropTypes.number,
  labels: PropTypes.array,
  placement: PropTypes.oneOf(['left', 'bottom']),
  labelHeight: PropTypes.number,
  labelWidth: PropTypes.number,
  labelPadding: PropTypes.number,
  overflowTop: PropTypes.number,
  overflowBottom: PropTypes.number,
  overflowLeft: PropTypes.number,
  overflowRight: PropTypes.number,
  tickLength: PropTypes.number,
  dtype: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

PlotAxis.defaultProps = {
  labelHeight: 12,
  tickLength: 6,
  labelPadding: 3,
  scale: 'linear',
  scientific: true, // Whether to use scientific notation, e.g. 1e+3
  siPostfix: false, // Whether to use SI postfixes, e.g. K, M, B
  decimals: 3 // How many decimals to show for custom labels
}

export default PlotAxis
