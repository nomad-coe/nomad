
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
import { makeStyles, useTheme, withStyles } from '@material-ui/core/styles'
import {
  Select,
  MenuItem,
  OutlinedInput
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import FilterChip from '../FilterChip'
import searchQuantities from '../../../searchQuantities'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import InputItem from './InputItem'
import { useSearchContext } from '../SearchContext'

// This forces the menu to have a fixed anchor instead of jumping around
const MenuProps = {
  getContentAnchorEl: null
}

// Customized input component
export const CustomInput = withStyles((theme) => ({
  input: {
    padding: theme.spacing(1),
    minHeight: '2.5rem'
  }
}))(OutlinedInput)

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  menuItem: {
    height: '2.1rem'
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
  count: {
    marginTop: theme.spacing(0.5),
    marginRight: theme.spacing(1),
    display: 'flex',
    justifyContent: 'flex-end',
    alignItems: 'center',
    color: theme.palette.text.disabled
  }
}))
const InputSelect = React.memo(({
  label,
  quantity,
  description,
  visible,
  initialScale,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const {useAgg, useFilterState, useFilterLocked} = useSearchContext()
  const agg = useAgg(quantity, visible)
  const [scale, setScale] = useState(initialScale)
  const [filter, setFilter] = useFilterState(quantity)
  const locked = useFilterLocked(quantity)
  const disabled = locked || (!(agg?.data && agg.data.length > 0))

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const title = label || def?.name

  const handleChange = useCallback((event) => {
    setFilter(new Set(event.target.value))
    event.preventDefault()
  }, [setFilter])

  // Create a memoized list of options
  const items = useMemo(() => {
    const items = []
    if (agg?.data) {
      for (let option of agg.data) {
        const value = option.value
        if (option.count > 0) {
          items.push(<MenuItem
            key={value}
            value={value}
            className={styles.menuItem}
            classes={{selected: styles.menuItemSelected}}
          >
            <InputItem
              key={value}
              value={value}
              selected={filter ? filter.has(value) : false}
              total={agg.total}
              variant="checkbox"
              count={option.count}
              scale={scale}
            />
          </MenuItem>)
        }
      }
    }
    return [items, agg?.data?.length || 0]
  }, [agg, filter, scale, styles])

  return <InputTooltip locked={locked} disabled={disabled}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <InputHeader
        quantity={quantity}
        label={title}
        description={desc}
        scale={scale}
        onChangeScale={setScale}
        disableAggSize
      />
      <Select
        disabled={disabled}
        multiple
        displayEmpty
        value={filter ? [...filter] : []}
        onChange={handleChange}
        input={<CustomInput/>}
        MenuProps={MenuProps}
        renderValue={(selected) => {
          return selected.length > 0
            ? <div className={styles.chips}>
              {selected.map((value) => <FilterChip
                locked={locked}
                key={value}
                label={value}
                color="primary"
              />)}
            </div>
            : <div className={styles.placeholder}>Click to select</div>
        }}
      >
        {items}
      </Select>
    </div>
  </InputTooltip>
})

InputSelect.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  visible: PropTypes.bool.isRequired,
  initialScale: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

InputSelect.defaultProps = {
  initialScale: 1
}

export default InputSelect
