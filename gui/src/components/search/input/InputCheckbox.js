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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  FormControlLabel,
  Checkbox,
  Tooltip,
  Typography
} from '@material-ui/core'
import PropTypes from 'prop-types'
import { isNil } from 'lodash'
import clsx from 'clsx'
import { useSearchContext, getValue } from '../SearchContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  }
}))

/**
 * A checkbox where the selection controls a true/false value within a quantity.
 */
const InputCheckbox = React.memo(({
  quantity,
  label,
  description,
  initialValue,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const {filterData, useFilterState, useFilterLocked} = useSearchContext()
  const [filter, setFilter] = useFilterState(quantity)
  const filterLocked = useFilterLocked(quantity)

  // Determine the description and units
  const def = filterData[quantity]
  const descFinal = description || def?.description || ''
  const labelFinal = label || def?.label
  const locked = !isNil(filterLocked) && def.global
  const val = getValue(def, filter, filterLocked, initialValue)

  const handleChange = useCallback((event, value) => {
    setFilter(value)
  }, [setFilter])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <Tooltip title={descFinal}>
      <FormControlLabel
        control={<Checkbox
          color="secondary"
          disabled={locked}
          checked={val}
          onChange={handleChange}
        />}
        label={<Typography>{labelFinal}</Typography>}
      />
    </Tooltip>
  </div>
})

InputCheckbox.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  initialValue: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default InputCheckbox

/**
 * A checkbox where the selection activates a specific value within a quantity.
 */
const useInputCheckboxValueStyles = makeStyles(theme => ({
  root: {
    '& > :last-child': {
      marginRight: -3
    }
  },
  control: {
    padding: 6
  }
}))
export const InputCheckboxValue = React.memo(({
  quantity,
  description,
  value,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useInputCheckboxValueStyles({classes: classes, theme: theme})
  const {filterData, useFilterState} = useSearchContext()
  const [filter, setFilter] = useFilterState(quantity)

  // Determine the description and units
  const def = filterData[quantity]
  const descFinal = description || def?.description || ''

  const handleChange = useCallback(() => {
    setFilter(old => {
      const newValue = old ? new Set(old) : new Set()
      newValue.has(value) ? newValue.delete(value) : newValue.add(value)
      return newValue
    })
  }, [setFilter, value])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <Tooltip title={descFinal}>
      <Checkbox
        color="secondary"
        checked={filter ? filter.has(value) : false}
        onChange={handleChange}
        className={styles.control}
        size="small"
      />
    </Tooltip>
  </div>
})

InputCheckboxValue.propTypes = {
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  value: PropTypes.any,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}
