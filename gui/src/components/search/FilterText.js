
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
import React, { useCallback, useMemo, useState } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { TextField, CircularProgress } from '@material-ui/core'
import Autocomplete, { createFilterOptions } from '@material-ui/lab/Autocomplete'
import parse from 'autosuggest-highlight/parse'
import match from 'autosuggest-highlight/match'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { Unit } from '../../units'
import { useApi } from '../apiV1'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
import { useFilterState } from './FilterContext'

const filterOptions = createFilterOptions({
  stringify: option => option,
  trim: true
})

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
  // const inputRef = useRef()
  const hasSuggestions = true
  const [suggestions, setSuggestions] = useState([])
  const [loading, setLoading] = useState(false)
  const api = useApi()

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const unitSI = def?.unit
  const unitLabel = useMemo(() => {
    const unit = unitSI && new Unit(unitSI)
    return unit && unit.label(units)
  }, [unitSI, units])
  const title = unitLabel ? `${name} (${unitLabel})` : name

  // Attach the filter hook
  const [filter, setFilter] = useFilterState(quantity, set)

  // Triggered when a value is picked by the user either by clicking a value or
  // pressing enter.
  const handleChange = useCallback((event, value, reason) => {
    value = value?.trim()
    setFilter(value)
  }, [setFilter])

  // When focus is lost on the element, the current input value is set as the
  // filter value. Notice that the 'autoSelect' property is not working to
  // achieve this: it will select the most recently highlighted value even if it
  // was not clicked.
  const handleClose = useCallback((event, reason) => {
    const value = event.target.value?.trim()
    setFilter(value)
  }, [setFilter])

  // Handle typing events. After a debounce time has expired, a list of
  // suggestion will be retrieved if they are available for this metainfo and
  // the input is deemed meaningful.
  const handleInputChange = useCallback((event, value, reason) => {
    if (!hasSuggestions) {
      return
    }
    value = value?.trim()
    if (!value || value.length < 2 || reason !== 'input') {
      setSuggestions([])
      return
    }
    setLoading(true)
    api.suggestions([quantity], value)
      .then(data => {
        setSuggestions(data)
      })
      .finally(() => setLoading(false))
  }, [quantity, api, hasSuggestions])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <FilterLabel label={title} description={desc}/>
    <Autocomplete
      freeSolo
      fullWidth
      value={filter === undefined ? null : filter}
      filterOptions={filterOptions}
      options={suggestions.map(option => option.value)}
      onInputChange={handleInputChange}
      onClose={handleClose}
      onChange={handleChange}
      getOptionSelected={(option, value) => false}
      renderOption={(option, { inputValue }) => {
        const matches = match(option, inputValue)
        const parts = parse(option, matches)
        return (
          <div>
            {parts.map((part, index) => (
              <span key={index} style={{ fontWeight: part.highlight ? 700 : 400 }}>
                {part.text}
              </span>
            ))}
          </div>
        )
      }}
      renderInput={(params) => (
        <TextField
          {...params}
          className={styles.textField}
          variant="outlined"
          InputProps={{
            ...params.InputProps,
            classes: {input: styles.input},
            endAdornment: (
              <React.Fragment>
                {loading ? <CircularProgress color="inherit" size={20} /> : null}
                {params.InputProps.endAdornment}
              </React.Fragment>
            )
          }}
        />
      )}
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
