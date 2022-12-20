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
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import InputHeader from './InputHeader'
import InputItem from './InputItem'
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
const InputRadio = React.memo(({
  quantity,
  label,
  description,
  initialValue,
  options,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const { filterData, useFilterState, useFilterLocked } = useSearchContext()
  const styles = useStyles({classes: classes, theme: theme})
  const [filter, setFilter] = useFilterState(quantity)
  const filterLocked = useFilterLocked(quantity)

  // Determine the description and units
  const def = filterData[quantity]
  const locked = !isNil(filterLocked) && filterData[quantity].global
  const descFinal = description || def?.description || ''
  const labelFinal = label || def?.label
  const val = getValue(def, filter, filterLocked, initialValue)

  const handleChange = useCallback((event, key, selected) => {
    setFilter(key)
  }, [setFilter])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <InputHeader
      quantity={quantity}
      label={labelFinal}
      description={descFinal}
      disableWidget
      disableStatistics
    />
    {options && Object.entries(options).map(([key, value]) =>
      <InputItem
        key={key}
        value={key}
        label={value.label}
        disabled={locked || value.disabled}
        selected={val === key}
        onChange={handleChange}
        tooltip={value.tooltip}
        variant="radio"
        disableStatistics
      />
    )}
  </div>
})

InputRadio.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  initialValue: PropTypes.string,
  options: PropTypes.object, // Mapping from option name to show label and tooltip
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default InputRadio
