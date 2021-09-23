
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
import {
  TextField,
  CircularProgress,
  Tooltip,
  IconButton
} from '@material-ui/core'
import Autocomplete from '@material-ui/lab/Autocomplete'
import CloseIcon from '@material-ui/icons/Close'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { Unit } from '../../../units'
import { useApi } from '../../api'
import searchQuantities from '../../../searchQuantities'
import InputLabel from './InputLabel'
import InputTooltip from './InputTooltip'
import { useSetFilter, useFilterLocked } from '../SearchContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  input: {
    padding: theme.spacing(1),
    height: '2.5rem'
  }
}))
const InputText = React.memo(({
  label,
  quantity,
  description,
  autocomplete,
  className,
  classes,
  units,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [suggestions, setSuggestions] = useState([])
  const [loading, setLoading] = useState(false)
  const [inputValue, setInputValue] = useState('')
  const [highlighted, setHighlighted] = useState({value: ''})
  const [open, setOpen] = useState(false)
  const [error, setError] = useState(false)
  const {api} = useApi()

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
  const setFilter = useSetFilter(quantity)
  const locked = useFilterLocked(quantity)

  const filterOptions = useCallback((options, {inputValue}) => {
    const trimmed = inputValue.trim().toLowerCase()
    return options.filter(option => {
      // If suggestionMode = suggestions, the results do not need to be filtered
      // (especially important because of fuzzy matches)
      if (autocomplete === 'suggestions') {
        return true
      }
      // Underscore can be replaced by a whitespace
      const optionClean = option.value.trim().toLowerCase()
      const matchUnderscore = optionClean.includes(trimmed)
      const matchNoUnderscore = optionClean.replace(/_/g, ' ').includes(trimmed)
      return matchUnderscore || matchNoUnderscore
    })
  }, [autocomplete])

  // Triggered when a value is submitted by pressing enter or clicking the
  // search icon.
  const handleSubmit = useCallback(() => {
    if (inputValue.trim().length === 0) {
      return
    }
    let valid = true
    if (valid) {
      // Submit to search context on successful validation.
      setFilter(inputValue)
      setInputValue('')
      setOpen(false)
    } else {
      setError(`Invalid query`)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [inputValue])

  const handleHighlight = useCallback((event, value, reason) => {
    setHighlighted(value)
  }, [])

  // Used to retrieve suggestions for this field.
  const fetchSuggestions = useCallback((input) => {
    input = input?.trim()
    if (input?.length > 0) {
      if (autocomplete === 'suggestions') {
        setLoading(true)
        api.suggestions([quantity], input)
          .then(data => {
            let res = []
            const esSuggestions = data[quantity]
            if (esSuggestions) {
              res = res.concat(esSuggestions.map(suggestion => ({
                value: suggestion.value
              })))
            }
            setSuggestions(res)
            setOpen(true)
          })
          .finally(() => setLoading(false))
      } else if (autocomplete === 'aggregations') {
      }
    }
  }, [api, autocomplete, quantity])

  // Submit the input value upon the element losing focus. This ensures that any
  // filters that are typed and not manually submitted (using enter) will still
  // either raise an error or be applied.
  const handleBlur = useCallback(() => {
    handleSubmit()
  }, [handleSubmit])

  // Handle clear button
  const handleClose = useCallback(() => {
    setInputValue('')
  }, [])

  // When enter is pressed, select currently highlighted value and close menu,
  // or if menu is not open submit the value.
  const handleEnter = useCallback((event) => {
    if (event.key === 'Enter') {
      if (open && highlighted?.value) {
        setInputValue(highlighted.value)
        setOpen(false)
      } else {
        handleSubmit()
      }
      event.stopPropagation()
      event.preventDefault()
    }
  }, [open, highlighted, handleSubmit])

  // Handle typing events. After a debounce time has expired, a list of
  // suggestion will be retrieved if they are available for this metainfo and
  // the input is deemed meaningful.
  const handleInputChange = useCallback((event, value, reason) => {
    setError(error => error ? undefined : null)
    setInputValue(value)

    // Suggestions are only retrieved in user input, not on clearing or
    // selecting a value
    if (value.trim() === '') {
      setSuggestions([])
      return
    }
    if (reason === 'input') {
      fetchSuggestions(value)
    }
  }, [fetchSuggestions])

  return <InputTooltip locked={locked}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <InputLabel
        quantity={quantity}
        label={title}
        description={desc}
        disableStatistics
      />
      <Autocomplete
        freeSolo
        disabled={locked}
        clearOnBlur={false}
        inputValue={inputValue}
        value={null}
        open={open}
        onOpen={() => setOpen(true)}
        onClose={() => setOpen(false)}
        fullWidth
        disableClearable
        classes={{endAdornment: styles.endAdornment}}
        filterOptions={filterOptions}
        options={suggestions}
        onBlur={handleBlur}
        onInputChange={handleInputChange}
        onHighlightChange={handleHighlight}
        getOptionLabel={option => option.value}
        getOptionSelected={(option, value) => false}
        renderInput={(params) => (
          <TextField
            {...params}
            variant="outlined"
            placeholder=""
            label={error || undefined}
            error={!!error}
            onKeyDown={handleEnter}
            InputLabelProps={{ shrink: true }}
            InputProps={{
              ...params.InputProps,
              endAdornment: (<>
                {loading ? <CircularProgress color="inherit" size={20} /> : null}
                {(inputValue?.length || null) && <>
                  <Tooltip title="Clear">
                    <IconButton
                      size="small"
                      onClick={handleClose}
                      className={styles.iconButton}
                      aria-label="clear"
                    >
                      <CloseIcon/>
                    </IconButton>
                  </Tooltip>
                </>}
              </>)
            }}
          />
        )}
      />
    </div>
  </InputTooltip>
})

InputText.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  description: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  autocomplete: PropTypes.string, // Determines the autocompletion mode: either 'aggregations', 'suggestions', or 'off'
  'data-testid': PropTypes.string
}

InputText.defaultProps = {
  autocomplete: 'suggestions'
}

export default InputText
