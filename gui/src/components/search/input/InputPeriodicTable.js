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
import React, { useMemo, useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import elementData from '../../../elementData'
import { useResizeDetector } from 'react-resize-detector'
import { Tooltip } from '@material-ui/core'
import InputHeader from './InputHeader'
import InputCheckbox from './InputCheckbox'
import { makeStyles, useTheme, lighten } from '@material-ui/core/styles'
import { useSearchContext } from '../SearchContext'
import { approxInteger } from '../../../utils'
import { getScaler } from '../../plotting/common'

/**
 * A single element in the periodic table.
*/
const useElementStyles = makeStyles(theme => ({
  text: {
    fontFamily: theme.typography.fontFamily,
    fill: theme.palette.text.default,
    userSelect: 'none',
    pointerEvents: 'none'
  },
  rect: {
    transition: 'fill .15s ease',
    cursor: 'pointer',
    '&:hover': {
      filter: 'drop-shadow(0px 0px 2px rgb(0 0 0 / 0.4))'
    }
  },
  rectDisabled: {
    cursor: 'auto',
    '&:hover': {
      filter: 'unset'
    }
  },
  label: {
    fontSize: '1rem',
    fontWeight: 600
  },
  number: {
    fontSize: '0.58rem'
  },
  count: {
    fontSize: '0.58rem'
  }
}))

const Element = React.memo(({
  x,
  y,
  width,
  height,
  element,
  selected,
  disabled,
  onClick,
  disableStatistics,
  max,
  count,
  scale
}) => {
  const styles = useElementStyles()
  const theme = useTheme()
  const scaler = useMemo(() => getScaler(scale, undefined, [0.2, 1]), [scale])
  const finalCount = useMemo(() => approxInteger(count || 0), [count])
  const finalScale = useMemo(() => scaler(count / max) || 0, [count, max, scaler])
  const disabledFinal = disabled && !selected
  const color = selected
    ? theme.palette.secondary.main
    : disabled
      ? '#eee'
      : disableStatistics
        ? theme.palette.primary.light
        : lighten(theme.palette.primary.light, 1 - finalScale)
  const strokeColor = disabledFinal
    ? theme.palette.text.disabled
    : '#555'
  const textColor = selected
    ? 'white'
    : disabled
      ? theme.palette.text.disabled
      : theme.palette.text.default

  const handleClick = useCallback(() => {
    if (disabledFinal) return
    onClick && onClick()
  }, [onClick, disabledFinal])

  return <Tooltip title={element.name}>
    <g>
      <rect
        data-testid={element.name}
        x={x}
        y={y}
        width={width}
        height={height}
        stroke={strokeColor}
        strokeWidth={1}
        fill={color}
        className={clsx(styles.rect, disabledFinal && styles.rectDisabled)}
        onClick={handleClick}
      />
      <text
        x={x + 0.5 * width}
        y={y + 0.53 * height}
        fill={textColor}
        textAnchor="middle"
        dominantBaseline="middle"
        className={clsx(styles.text, styles.label)}>{element.symbol}
      </text>
      <text
        x={x + 0.04 * width}
        y={y + 0.03 * height}
        fill={textColor}
        textAnchor="start"
        dominantBaseline="hanging"
        className={clsx(styles.text, styles.number)}>{element.number}
      </text>
      {!disableStatistics && <text
        x={x + 0.95 * width}
        y={y + 0.93 * height}
        fill={textColor}
        textAnchor="end"
        dominantBaseline="auto"
        className={clsx(styles.text, styles.number)}>{finalCount}
      </text>}
    </g>
  </Tooltip>
})

Element.propTypes = {
  x: PropTypes.number.isRequired,
  y: PropTypes.number.isRequired,
  width: PropTypes.number.isRequired,
  height: PropTypes.number.isRequired,
  element: PropTypes.object.isRequired,
  max: PropTypes.number,
  count: PropTypes.number,
  scale: PropTypes.string,
  selected: PropTypes.bool,
  disabled: PropTypes.bool,
  onClick: PropTypes.func,
  disableStatistics: PropTypes.bool
}

/**
 * Displays an interactive periodic table.
*/
const usePeriodicTableStyles = makeStyles(theme => ({
  root: {
    minHeight: 0, // added min-height: 0 to prevent relayouting when within flexbox
    flexGrow: 1,
    width: '100%',
    height: '100%',
    position: 'relative'
  },
  container: {
    position: 'absolute',
    top: theme.spacing(-0.2),
    left: '10%',
    textAlign: 'center'
  }
}))
export const PeriodicTable = React.memo(({
  quantity,
  visible,
  scale,
  anchored,
  disableStatistics,
  aggId,
  className,
  'data-testid': testID
}) => {
  const styles = usePeriodicTableStyles()
  const { height, width, ref } = useResizeDetector()
  const {useFilterState, useAgg, useIsStatisticsEnabled} = useSearchContext()
  const isStatisticsEnabled = useIsStatisticsEnabled()
  const [filter, setFilter] = useFilterState(quantity)
  const aggConfig = useMemo(() => ({type: 'terms'}), [])
  const agg = useAgg(quantity, visible, aggId, aggConfig)
  const availableValues = useMemo(() => {
    const elementCountMap = {}
    agg?.data && agg.data.forEach((value) => { elementCountMap[value.value] = value.count })
    return elementCountMap
  }, [agg])
  disableStatistics = anchored
    ? false
    : isNil(disableStatistics) ? !isStatisticsEnabled : disableStatistics

  const onElementClicked = useCallback((element) => {
    setFilter(old => {
      let newValues
      if (old) {
        const isSelected = old?.has(element)
        isSelected ? old.delete(element) : old.add(element)
        newValues = new Set(old)
      } else {
        newValues = new Set([element])
      }
      return newValues
    })
  }, [setFilter])

  const max = useMemo(() => {
    return agg
      ? Math.max(...agg.data
        .filter(option => !filter?.has(option.value))
        .map(option => option.count))
      : 0
  }, [agg, filter])

  const padding = 1
  const margin = 3
  const cellWidth = (width - 2 * padding - 17 * margin) / 18
  const cellHeight = (height - 2 * padding - 9 * margin) / 10

  return <div ref={ref} className={clsx(styles.root, className)} data-testid={testID}>
    <svg
      width="100%"
      height="100%"
      viewBox={`0 0 ${width || 0} ${height || 0}`}
      shapeRendering='optimizeSpeed'
      textRendering='optimizeSpeed'
    >
      {(width && height) &&
        elementData.elements.map((element) => (
          <Element
            key={element.name}
            element={element}
            disabled={!availableValues[element.symbol]}
            onClick={() => onElementClicked(element.symbol)}
            selected={filter?.has(element.symbol)}
            max={max}
            count={availableValues[element.symbol]}
            scale={scale}
            x={padding + (element.xpos - 1) * (cellWidth + margin)}
            y={padding + (element.ypos - 1) * (cellHeight + margin)}
            width={cellWidth}
            height={cellHeight}
            disableStatistics={disableStatistics}
          />
        ))
      }
    </svg>
    <div className={styles.container}>
      <InputCheckbox
        quantity="exclusive"
        label="only compositions that exclusively contain these atoms"
      ></InputCheckbox>
    </div>
  </div>
})

PeriodicTable.propTypes = {
  quantity: PropTypes.string,
  label: PropTypes.string,
  description: PropTypes.string,
  visible: PropTypes.bool,
  scale: PropTypes.string,
  anchored: PropTypes.bool,
  disableStatistics: PropTypes.bool,
  aggId: PropTypes.string,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

PeriodicTable.defaultProps = {
  aggId: 'default',
  'data-testid': 'periodictable'
}

/**
 * Periodic table that directly shows the targeted quantity along with some
 * extra controls.
 */
const useInputPeriodicTableStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'column'
  }
}))
const InputPeriodicTable = React.memo(({
    quantity,
    label,
    description,
    visible,
    initialScale,
    anchored,
    disableStatistics,
    aggId,
    className
  }) => {
  const styles = useInputPeriodicTableStyles()
  const [scale, setScale] = useState(initialScale)
  const {filterData} = useSearchContext()

  // Determine the description and title
  const def = filterData[quantity]
  const descFinal = description || def?.description || ''
  const labelFinal = label || def?.label

  return <div className={clsx(styles.root, className)}>
    <InputHeader
      quantity={quantity}
      label={labelFinal}
      description={descFinal}
      scale={scale}
      onChangeScale={setScale}
      disableStatistics={disableStatistics}
      disableAggSize
      anchored={anchored}
    />
    <PeriodicTable
      quantity={quantity}
      visible={visible}
      scale={scale}
      anchored={anchored}
      disableStatistics={disableStatistics}
      aggId={aggId}
    />
  </div>
})

InputPeriodicTable.propTypes = {
  quantity: PropTypes.string,
  label: PropTypes.string,
  description: PropTypes.string,
  visible: PropTypes.bool,
  initialScale: PropTypes.string,
  anchored: PropTypes.bool,
  disableStatistics: PropTypes.bool,
  aggId: PropTypes.string,
  className: PropTypes.string
}

InputPeriodicTable.defaultProps = {
  initialScale: 'linear'
}

export default InputPeriodicTable
