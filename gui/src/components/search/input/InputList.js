
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
import React, { useState, useCallback, useMemo } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Typography } from '@material-ui/core'
import PropTypes from 'prop-types'
import { useResizeDetector } from 'react-resize-detector'
import clsx from 'clsx'
import searchQuantities from '../../../searchQuantities'
import InputLabel from './InputLabel'
import InputTooltip from './InputTooltip'
import InputItem, { inputItemHeight } from './InputItem'
import { useFilterState, useFilterLocked, useAgg } from '../SearchContext'

/**
 * Displays a list of options with fixed maximum size. Only options that are
 * present in the current search results are displayed. The options are sorted
 * by occurrence and the number of displayed results can be changed.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  menuItem: {
    height: inputItemHeight
  },
  menuItemSelected: {
    '&.Mui-selected': {
      backgroundColor: 'transparent'
    }
  },
  chips: {
    display: 'flex',
    flexWrap: 'wrap'
  },
  placeholder: {
    marginLeft: theme.spacing(0.5),
    height: '2.5rem',
    display: 'flex',
    alignItems: 'center',
    color: theme.palette.text.disabled
  },
  icon: {
    right: theme.spacing(1)
  },
  spacer: {
    overflow: 'hidden',
    flex: 1
  },
  count: {
    marginTop: theme.spacing(0.5),
    marginRight: theme.spacing(1),
    display: 'flex',
    justifyContent: 'flex-end',
    alignItems: 'center',
    color: theme.palette.text.disabled
  }
}))
const InputList = React.memo(({
  label,
  quantity,
  description,
  visible,
  initialScale,
  draggable,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [aggSize, setAggSize] = useState(0)
  const agg = useAgg(quantity, visible, false)
  const [scale, setScale] = useState(initialScale)
  const [filter, setFilter] = useFilterState(quantity)
  const locked = useFilterLocked(quantity)
  const disabled = locked || (!(agg?.data && agg.data.length > 0))
  const { height, ref } = useResizeDetector()

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const title = label || def?.name

  const handleChange = useCallback((key, value) => {
    setFilter(old => {
      if (!old) return new Set([key])
      const newValue = new Set(old)
      value ? newValue.add(key) : newValue.delete(key)
      return newValue
    })
  }, [setFilter])

  // Create a memoized list of options
  const [items, nItems, nTotal] = useMemo(() => {
    const items = []
    let index = 0
    let nItems = Math.floor(height / inputItemHeight)
    if (agg?.data) {
      for (let option of agg.data) {
        const value = option.value
        if (option.count > 0) {
          if (index < nItems) {
            items.push(<InputItem
              key={value}
              value={value}
              selected={filter ? filter.has(value) : false}
              total={agg.total}
              onChange={handleChange}
              variant="checkbox"
              count={option.count}
              scale={scale}
            />)
          }
          ++index
        }
      }
    }
    return [items, nItems, index]
  }, [agg, filter, scale, handleChange, height])

  return <InputTooltip locked={locked} disabled={disabled}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <InputLabel
        quantity={quantity}
        label={title}
        description={desc}
        scale={scale}
        onChangeScale={setScale}
        aggSize={aggSize}
        onChangeAggSize={setAggSize}
        draggable={draggable}
      />
      <div ref={ref} className={styles.spacer}>
        {items}
      </div>
      <div className={styles.count}>
        <Typography variant="overline">{`${Math.min(nItems, nTotal)}/${nTotal}`}</Typography>
      </div>
    </div>
  </InputTooltip>
})

InputList.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  visible: PropTypes.bool.isRequired,
  initialScale: PropTypes.number,
  draggable: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

InputList.defaultProps = {
  initialScale: 1
}

export default InputList
