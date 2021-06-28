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
import React, { useState, useCallback, useEffect } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Tooltip
} from '@material-ui/core'
import {
  KeyboardDatePicker
} from '@material-ui/pickers'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
import { useAgg, useFilterState } from './FilterContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  row: {
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    justifyContent: 'flex-start'
  },
  dash: {
    padding: theme.spacing(1)
  },
  date: {
    flex: 1
  }
}))
const FilterDate = React.memo(({
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
  const [startDate, setStartDate] = useState(new Date())
  const [endDate, setEndDate] = useState(new Date())
  // eslint-disable-next-line no-unused-vars
  const [filter, setFilter] = useFilterState(quantity)
  const [startGlobal, endGlobal] = useAgg(quantity, 'min_max', true, visible)
  const disabled = false

  // The range is automatically adjusted according to global min/max of the
  // field
  useEffect(() => {
    if (!isNil(endGlobal) && !isNil(startGlobal)) {
      console.log(startGlobal, endGlobal)
      console.log(typeof startGlobal)
      if (filter === undefined) {
        setStartDate(new Date(startGlobal))
        setEndDate(new Date(endGlobal))
      }
    }
  }, [startGlobal, endGlobal, filter])

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const title = label || def?.name

  const handleDateChange = useCallback((date) => {
  }, [])

  return <Tooltip title={disabled ? 'No values available with current query.' : ''}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <FilterLabel label={title} description={desc}/>
      <div className={styles.row}>
        <KeyboardDatePicker
          className={styles.date}
          autoOk
          disabled={disabled}
          variant="inline"
          inputVariant="outlined"
          label="Start date"
          format="dd/MM/yyyy"
          value={startDate}
          InputAdornmentProps={{ position: 'start' }}
          onChange={date => handleDateChange(date)}
        />
        <div className={styles.dash}>-</div>
        <KeyboardDatePicker
          className={styles.date}
          autoOk
          disabled={disabled}
          variant="inline"
          inputVariant="outlined"
          label="End date"
          format="dd/MM/yyyy"
          value={endDate}
          InputAdornmentProps={{ position: 'start' }}
          onChange={date => handleDateChange(date)}
        />
      </div>
    </div>
  </Tooltip>
})

FilterDate.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  visible: PropTypes.bool.isRequired,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterDate
