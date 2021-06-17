
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
import React, { useCallback } from 'react'
import { makeStyles, useTheme, withStyles } from '@material-ui/core/styles'
import {
  Select,
  MenuItem,
  Checkbox,
  ListItemText,
  OutlinedInput,
  Tooltip
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
import FilterChip from './FilterChip'
import { useFilterState, useAgg } from './FilterContext'

// This forces the menu to have a fixed anchor instead of jumping around
const MenuProps = {
  getContentAnchorEl: null
}

// Customized input component
const CustomInput = withStyles((theme) => ({
  input: {
    padding: theme.spacing(1),
    minHeight: '2.5rem'
  }
}))(OutlinedInput)

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  select: {
    width: '100%'
  },
  chips: {
    display: 'flex',
    flexWrap: 'wrap'
  },
  icon: {
    right: theme.spacing(1)
  }
}))
const FilterSelect = React.memo(({
  label,
  quantity,
  description,
  visible,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const options = useAgg(quantity, 'terms', true, visible)
  const [filter, setFilter] = useFilterState(quantity)
  const disabled = !(options && options.length > 0)

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const title = label || def?.name

  // Create a list of options
  const menuItems = options && options
    .filter((option) => option.count !== 0)
    .map((option) => {
      const value = option.value
      return <MenuItem key={value} value={value}>
        <Checkbox checked={filter ? filter.has(value) : false} />
        <ListItemText primary={value} />
      </MenuItem>
    })

  const handleChange = useCallback((event) => {
    setFilter(new Set(event.target.value))
  }, [setFilter])

  return <Tooltip title={disabled ? 'No values available with current query.' : ''}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <FilterLabel label={title} description={desc}/>
      <Select
        disabled={disabled}
        multiple
        value={filter ? [...filter] : []}
        onChange={handleChange}
        className={styles.select}
        classes={{icon: styles.icon}}
        input={<CustomInput/>}
        MenuProps={MenuProps}
        renderValue={(selected) => (
          <div className={styles.chips}>
            {selected.map((value) => <FilterChip key={value} label={value}/>)}
          </div>
        )}
      >
        {menuItems}
      </Select>
    </div>
  </Tooltip>
})

FilterSelect.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  visible: PropTypes.bool.isRequired,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterSelect
