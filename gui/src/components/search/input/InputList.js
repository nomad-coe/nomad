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
import { useRecoilValue } from 'recoil'
import { isNil } from 'lodash'
import { useResizeDetector } from 'react-resize-detector'
import clsx from 'clsx'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import InputItem, { inputItemHeight } from './InputItem'
import InputUnavailable from './InputUnavailable'
import Placeholder from '../../visualization/Placeholder'
import { guiState } from '../../GUIMenu'
import { useSearchContext } from '../SearchContext'
import { pluralize } from '../../../utils'

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
  icon: {
    right: theme.spacing(1)
  },
  spacer: {
    flex: 1,
    minHeight: 0
  },
  count: {
    marginTop: theme.spacing(0.5),
    marginRight: theme.spacing(1),
    display: 'flex',
    justifyContent: 'flex-end',
    alignItems: 'center',
    color: theme.palette.text.disabled
  },
  placeholder: {
    top: theme.spacing(0),
    left: theme.spacing(0),
    right: theme.spacing(0),
    bottom: theme.spacing(0)
  }
}))
const InputList = React.memo(({
  label,
  quantity,
  description,
  visible,
  initialScale,
  draggable,
  aggId,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const {filterData, useAgg, useFilterState, useFilterLocked} = useSearchContext()
  const styles = useStyles({classes: classes, theme: theme})
  const [scale, setScale] = useState(initialScale || filterData[quantity].scale)
  const aggIndicator = useRecoilValue(guiState('aggIndicator'))
  const [filter, setFilter] = useFilterState(quantity)
  const locked = useFilterLocked(quantity)
  const { height, ref } = useResizeDetector()
  const aggSize = useMemo(() => Math.floor(height / inputItemHeight), [height])
  const agg = useAgg(quantity, !isNil(height) && visible, aggSize, aggId)
  const max = agg ? Math.max(...agg.data.map(option => option.count)) : 0

  // Determine the description and units
  const def = filterData[quantity]
  const descFinal = description || def?.description || ''
  const labelFinal = label || def?.label

  const handleChange = useCallback((event, key, selected) => {
    setFilter(old => {
      if (!old) return new Set([key])
      const newValue = new Set(old)
      selected ? newValue.add(key) : newValue.delete(key)
      return newValue
    })
  }, [setFilter])

  const [aggComp, nShown] = useMemo(() => {
    let aggComp
    let nShown = 0
    if (!agg && aggIndicator === 'on') {
      aggComp = <Placeholder
        variant="rect"
        classes={{placeholder: styles.placeholder}}
        data-testid={`${testID}-placeholder`}
      />
      nShown = 0
    } else {
      if (agg?.data && agg.data.length > 0) {
        aggComp = []
        const maxSize = Math.min(aggSize, agg.data.length)
        for (let i = 0; i < maxSize; ++i) {
          const option = agg.data[i]
          if (option.count > 0 && nShown < maxSize) {
            aggComp.push(<InputItem
              key={option.value}
              value={option.value}
              selected={filter ? filter.has(option.value) : false}
              max={max}
              onChange={handleChange}
              variant="checkbox"
              count={option.count}
              scale={scale}
              disabled={locked}
            />)
            ++nShown
          }
        }
      } else {
        aggComp = <InputUnavailable/>
        nShown = 0
      }
    }
    return [aggComp, nShown]
  }, [agg, aggIndicator, aggSize, filter, handleChange, locked, max, scale, styles, testID])

  return <InputTooltip locked={locked}>
    <div className={clsx(className, styles.root)}>
      <InputHeader
        quantity={quantity}
        label={labelFinal}
        description={descFinal}
        scale={scale}
        onChangeScale={setScale}
        draggable={draggable}
        locked={locked}
      />
      <div ref={ref} className={styles.spacer}>
        {aggComp}
      </div>
      {nShown !== 0 && <div className={styles.count}>
        <Typography variant="overline">
          {nShown === 1
            ? 'Showing the only item'
            : agg.data.length <= aggSize
              ? `Showing all ${pluralize('item', nShown, true)}`
              : `Showing top ${pluralize('item', nShown, true)}`
          }
        </Typography>
      </div>}
    </div>
  </InputTooltip>
})

InputList.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  visible: PropTypes.bool.isRequired,
  initialScale: PropTypes.string,
  draggable: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  aggId: PropTypes.string,
  'data-testid': PropTypes.string
}

InputList.defaultProps = {
  aggId: 'default'
}

export default InputList
