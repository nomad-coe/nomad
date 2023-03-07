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
import React, { useMemo, useCallback, useState, useEffect, useRef } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import elementData from '../../../elementData'
import {
  Typography,
  ButtonBase,
  Tooltip
} from '@material-ui/core'
import InputHeader from './InputHeader'
import InputCheckbox from './InputCheckbox'
import { makeStyles } from '@material-ui/core/styles'
import { useSearchContext } from '../SearchContext'
import { approxInteger } from '../../../utils'
import { getScaler } from '../../plotting/common'

// A fixed 2D, 10x18 array for the element data.
const elements = []
for (let i = 0; i < 10; i++) {
  elements[i] = Array.apply(null, Array(18))
}
elementData.elements.forEach(element => {
  elements[element.ypos - 1][element.xpos - 1] = element
  element.category = element.category.replace(' ', '')
})

/**
 * A single element in the periodic table.
*/
const useElementStyles = makeStyles(theme => ({
  root: {
    top: 1,
    bottom: 1,
    left: 1,
    right: 1,
    position: 'absolute'
  },
  fit: {
    top: 0,
    bottom: 0,
    left: 0,
    right: 0,
    position: 'absolute'
  },
  bg: {
    opacity: 0,
    willChange: 'opacity',
    transition: 'opacity 250ms',
    backgroundColor: theme.palette.primary.light
  },
  disabled: {
    opacity: 1,
    willChange: 'opacity',
    transition: 'opacity 250ms',
    backgroundColor: '#eee'
  },
  selected: {
    backgroundColor: theme.palette.secondary.main,
    display: 'none'
  },
  visible: {
    display: 'block'
  },
  button: {
    color: theme.palette.text.default,
    width: '100%',
    height: '100%',
    border: '1px solid',
    borderColor: '#555',
    textAlign: 'center',
    fontFamily: theme.typography.fontFamily,
    fontSize: '1rem',
    fontWeight: 600,
    '&:hover': {
      boxShadow: theme.shadows[4]
    }
  },
  buttonSelected: {
    color: 'white'
  },
  buttonDisabled: {
    borderColor: '#999',
    color: theme.palette.text.disabled
  },
  number: {
    position: 'absolute',
    top: 0,
    left: 2,
    margin: 0,
    padding: 0,
    fontSize: 9,
    pointerEvents: 'none'
  },
  count: {
    position: 'absolute',
    bottom: 0,
    right: 2,
    margin: 0,
    padding: 0,
    fontSize: 9,
    pointerEvents: 'none'
  },
  textSelected: {
    color: 'white'
  },
  textDisabled: {
    color: '#BDBDBD'
  },
  symbol: {
    marginTop: -2
  }
}))

const Element = React.memo(({
  element,
  selected,
  disabled,
  onClick,
  disableStatistics,
  max,
  count,
  scale,
  localFilter
}) => {
  const styles = useElementStyles()

  // Calculate the approximated count and the final scaled value
  const scaler = useMemo(() => getScaler(scale, undefined, [0.2, 1]), [scale])
  const finalCount = useMemo(() => approxInteger(count || 0), [count])
  const finalScale = useMemo(() => scaler(count / max) || 0, [count, max, scaler])

  // Dynamically calculated styles. The background color is formed by animating
  // opacity: opacity animation can be GPU-accelerated by the browser unlike
  // animating the color property.
  const useDynamicStyles = makeStyles((theme) => {
    return {
      bg: { opacity: disableStatistics
        ? 0.4
        : (isNil(count) || isNil(max))
          ? 0
          : finalScale
      },
      disabled: { opacity: disabled ? 1 : 0 }
    }
  })
  const dynamicStyles = useDynamicStyles()

  const [selectedInternal, setSelectedInternal] = useState(selected)
  useEffect(() => {
    setSelectedInternal(selected)
  }, [selected])
  const disabledInternal = selectedInternal ? false : disabled

  const handleClick = useCallback(() => {
    setSelectedInternal(old => {
      const newValue = !old
      if (newValue) {
        localFilter.add(element.symbol)
      } else {
        localFilter.delete(element.symbol)
      }
      return newValue
    })
    onClick()
  }, [onClick, element, localFilter])

  return <div className={styles.root}>
    <div className={clsx(styles.fit, styles.bg, dynamicStyles.bg)}/>
    <div className={clsx(styles.fit, styles.disabled, dynamicStyles.disabled)}/>
    <div className={clsx(styles.fit, styles.selected, selectedInternal && styles.visible)}/>
    <Tooltip title={element.name}>
      <span className={styles.fit}>
        <ButtonBase
          className={clsx(
            styles.fit,
            styles.button,
            selectedInternal && styles.buttonSelected,
            disabledInternal && styles.buttonDisabled)
          }
          disabled={disabledInternal}
          onClick={handleClick}
          variant="contained"
        >
          <span className={styles.symbol}>{element.symbol}</span>
        </ButtonBase>
      </span>
    </Tooltip>
    <Typography
      className={clsx(
        styles.number,
        selectedInternal && styles.textSelected,
        disabledInternal && styles.textDisabled
      )}
      variant="caption"
    >
      {element.number}
    </Typography>
    {(!disableStatistics) && <Typography
      className={clsx(
        styles.count,
        selectedInternal && styles.textSelected,
        disabledInternal && styles.textDisabled
      )}
      variant="caption"
    >
      {finalCount}
    </Typography>}
  </div>
})

Element.propTypes = {
  element: PropTypes.object.isRequired,
  onClick: PropTypes.func,
  selected: PropTypes.bool,
  disabled: PropTypes.bool,
  disableStatistics: PropTypes.bool,
  max: PropTypes.number,
  count: PropTypes.number,
  scale: PropTypes.string,
  localFilter: PropTypes.object
}

/**
 * Represents a single element in the periodic table.
*/
const useTableStyles = makeStyles(theme => ({
  root: {
    height: '100%',
    width: '100%',
    display: 'flex',
    flexDirection: 'column'
  },
  container: {
    flex: '1 1 100%',
    position: 'relative'
  },
  table: {
    width: '100%',
    height: '100%',
    borderSpacing: 0,
    tableLayout: 'fixed'
  },
  td: {
    position: 'relative'
  },
  formContainer: {
    position: 'absolute',
    top: theme.spacing(-0.2),
    left: '10%',
    textAlign: 'center'
  }
}))

function eqSet(as, bs) {
  if (isNil(as) || isNil(bs)) return false
  if (as.size !== bs.size) return false
  for (const a of as) if (!bs.has(a)) return false
  return true
}

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
  const styles = useTableStyles()
  const {useFilterState, useAgg, useIsStatisticsEnabled} = useSearchContext()
  const isStatisticsEnabled = useIsStatisticsEnabled()
  const [filter, setFilter] = useFilterState(quantity)
  const localFilter = useRef(new Set())
  const [update, setUpdate] = useState(0)
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

  // The selected state of the periodic filter is kept in a local reference.
  // This way simply selecting an element does not cause a full re-render of the
  // table. To handle external changes to the filter state, the local state is
  // synced each time a change is triggered and only if the states differ, a
  // re-render is issued.
  useEffect(() => {
    if (!eqSet(filter, localFilter.current)) {
      localFilter.current = new Set(filter)
      setUpdate(old => old + 1)
    }
  }, [filter, setUpdate, localFilter])

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

  const table = useMemo(() => {
    // The colors are normalized with respect to the maximum aggregation size
    // for an unselected element.
    const max = agg
      ? Math.max(...agg.data
        .filter(option => !localFilter.current.has(option.value))
        .map(option => option.count))
      : 0
    return <div className={clsx(styles.root, className)} data-testid={testID}>
      <div className={styles.container}>
        <table className={styles.table}>
          <tbody>
            {elements.map((row, i) => (
              <tr key={i}>
                {row.map((element, j) => (
                  <td key={j} className={styles.td}>
                    {element
                      ? <Element
                        element={element}
                        disabled={!availableValues[element.symbol]}
                        onClick={() => onElementClicked(element.symbol)}
                        selected={localFilter.current.has(element.symbol)}
                        max={max}
                        count={availableValues[element.symbol]}
                        localFilter={localFilter.current}
                        scale={scale}
                        disableStatistics={disableStatistics}
                      />
                      : null}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
        <div className={styles.formContainer}>
          <InputCheckbox
            quantity="exclusive"
            label="only compositions that exclusively contain these atoms"
          ></InputCheckbox>
        </div>
      </div>
    </div>
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [agg, availableValues, onElementClicked, styles, update, scale, disableStatistics, className, testID])

  return table
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
