
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
import React, { useRef, useCallback, useMemo } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { TextField } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { Unit } from '../../units'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
import { useSetFilter } from './FilterContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  textField: {
    marginTop: theme.spacing(1)
  },
  input: {
    padding: '16px 12px'
  }
}))
const FilterText = React.memo(({
  label,
  quantity,
  description,
  className,
  classes,
  units,
  set,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const inputRef = useRef()

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const unitSI = def?.unit
  const unit = useMemo(() => {
    return unitSI && new Unit(unitSI, units)
  }, [unitSI, units])
  const unitLabel = unit && unit.label()
  const title = unitLabel ? `${name} (${unitLabel})` : name

  // Attach the filter hook
  const setFilter = useSetFilter(quantity, set)

  // Handle input submission
  const handleKeyUp = useCallback(event => {
    if (event.key === 'Enter') {
      const value = inputRef.current.value.trim()
      if (value) {
        setFilter(inputRef.current.value)
        inputRef.current.value = ''
      }
    }
  }, [setFilter])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <FilterLabel label={title} description={desc}/>
    <TextField
      variant="outlined"
      fullWidth
      inputRef={inputRef}
      onKeyUp={handleKeyUp}
      className={styles.textField}
      InputProps={{classes: {input: styles.input}}}
    />
  </div>
})

FilterText.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  description: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string,
  set: PropTypes.object
}

export default FilterText
