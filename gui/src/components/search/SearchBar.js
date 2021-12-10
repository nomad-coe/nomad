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
import React, { useCallback, useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import Autocomplete from '@material-ui/lab/Autocomplete'
import { makeStyles } from '@material-ui/core/styles'
import SearchIcon from '@material-ui/icons/Search'
import CloseIcon from '@material-ui/icons/Close'
import {
  TextField,
  CircularProgress,
  Paper,
  Tooltip
} from '@material-ui/core'
import IconButton from '@material-ui/core/IconButton'
import { useUnits } from '../../units'
import { DType, getDatatype } from '../../utils'
import { useSuggestions } from '../../hooks'
import { useSearchContext, toGUIFilterSingle } from './SearchContext'
import searchQuantities from '../../searchQuantities'
import {
  filterFullnames,
  filterAbbreviations
} from './FilterRegistry'

const opMap = {
  '<=': 'lte',
  '>=': 'gte',
  '>': 'gt',
  '<': 'lt'
}
const opMapReverse = {
  '<=': 'gte',
  '>=': 'lte',
  '>': 'lt',
  '<': 'gt'
}

// Decides which options are shown
const filterOptions = (options, {inputValue}) => {
  const trimmed = inputValue.trim().toLowerCase()
  return options.filter(option => {
    // ES results do not need to be filtered at all
    const category = option.category
    if (category !== 'quantity name') {
      return true
    }
    // Underscore can be replaced by a whitespace
    const optionClean = option.value.trim().toLowerCase()
    const matchUnderscore = optionClean.includes(trimmed)
    const matchNoUnderscore = optionClean.replace(/_/g, ' ').includes(trimmed)
    return matchUnderscore || matchNoUnderscore
  })
}

// Customized paper component for the autocompletion options
const CustomPaper = (props) => {
  return <Paper elevation={3} {...props} />
}

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    position: 'relative'
  },
  notchedOutline: {
    borderColor: 'rgba(0, 0, 0, 0.0)'
  },
  iconButton: {
    padding: 10
  },
  divider: {
    height: '2rem'
  },
  endAdornment: {
    position: 'static'
  },
  examples: {
    position: 'absolute',
    left: 0,
    right: 0,
    top: 'calc(100% + 4px)',
    padding: theme.spacing(2),
    fontStyle: 'italic'
  }
}))

/**
 * This component shows a searchbar with autocomplete functionality. It does its
 * on API calls to provide autocomplete suggestion options.
 */
const SearchBar = React.memo(({
  className
}) => {
  const styles = useStyles()
  const units = useUnits()
  const {
    filters,
    filterData,
    useSetFilters,
    useFiltersLocked
  } = useSearchContext()
  const [inputValue, setInputValue] = useState('')
  const [suggestionInput, setSuggestionInput] = useState('')
  const [highlighted, setHighlighted] = useState({value: ''})
  const [open, setOpen] = useState(false)
  const [error, setError] = useState(false)
  const filtersLocked = useFiltersLocked()
  const setFilter = useSetFilters()
  const quantitySetSuggestable = useMemo(
    () => new Set([...filters].filter(q => searchQuantities[q]?.suggestion)),
    [filters]
  )
  const [quantityList, setQuantityList] = useState([...quantitySetSuggestable])
  const [suggestions, loading] = useSuggestions(quantityList, suggestionInput)
  const quantitySuggestions = useMemo(() => {
    const suggestions = []
    for (let q of filters) {
      if (!filterData[q].nested) {
        suggestions.push({
          value: filterAbbreviations[q] || q,
          category: 'quantity name'
        })
      }
    }
    return suggestions
  }, [filterData, filters])

  // Triggered when a value is submitted by pressing enter or clicking the
  // search icon.
  const handleSubmit = useCallback(() => {
    if (inputValue.trim().length === 0) {
      return
    }
    const reString = '[^\\s=<>](?:[^=<>]*[^\\s=<>])?'
    const op = '(?:<|>)=?'
    let valid = false
    let quantityFullname
    let queryValue

    // Equality query
    const equals = inputValue.match(new RegExp(`^\\s*(${reString})\\s*=\\s*(${reString})\\s*$`))
    if (equals) {
      const quantityName = equals[1]
      quantityFullname = filterFullnames[quantityName] || quantityName
      if (!filters.has(quantityFullname)) {
        setError(`Unknown quantity name`)
        return
      }
      try {
        queryValue = toGUIFilterSingle(quantityFullname, equals[2], units)
      } catch (error) {
        setError(`Invalid value for this metainfo. Please check your syntax.`)
        return
      }
      valid = true
    }

    // Simple LTE/GTE query
    if (!valid) {
      const ltegte = inputValue.match(new RegExp(`^\\s*(${reString})\\s*(${op})\\s*(${reString})\\s*$`))
      if (ltegte) {
        const a = ltegte[1]
        const op = ltegte[2]
        const b = ltegte[3]
        const aFullname = filterFullnames[a]
        const bFullname = filterFullnames[b]
        const isAQuantity = filters.has(aFullname)
        const isBQuantity = filters.has(bFullname)
        if (!isAQuantity && !isBQuantity) {
          setError(`Unknown quantity name`)
          return
        }
        quantityFullname = isAQuantity ? aFullname : bFullname
        const dtype = getDatatype(quantityFullname)
        if (dtype !== DType.Number && dtype !== DType.Timestamp) {
          setError(`Cannot perform range query for a non-numeric quantity.`)
          return
        }
        let quantityValue
        try {
          quantityValue = toGUIFilterSingle(quantityFullname, isAQuantity ? b : a, units)
        } catch (error) {
          console.log(error)
          setError(`Invalid value for this metainfo. Please check your syntax.`)
          return
        }
        queryValue = {}
        queryValue[opMap[op]] = quantityValue
        valid = true
      }
    }

    // Sandwiched LTE/GTE query
    if (!valid) {
      const ltegteSandwich = inputValue.match(new RegExp(`^\\s*(${reString})\\s*(${op})\\s*(${reString})\\s*(${op})\\s*(${reString})\\s*$`))
      if (ltegteSandwich) {
        const a = ltegteSandwich[1]
        const op1 = ltegteSandwich[2]
        const b = ltegteSandwich[3]
        const op2 = ltegteSandwich[4]
        const c = ltegteSandwich[5]
        quantityFullname = filterFullnames[b]
        const dtype = getDatatype(quantityFullname)
        if (dtype !== DType.Number && dtype !== DType.Timestamp) {
          setError(`Cannot perform range query for a non-numeric quantity.`)
          return
        }
        const isBQuantity = filters.has(quantityFullname)
        if (!isBQuantity) {
          setError(`Unknown quantity name`)
          return
        }

        queryValue = {}
        try {
          queryValue[opMapReverse[op1]] = toGUIFilterSingle(quantityFullname, a, units)
          queryValue[opMap[op2]] = toGUIFilterSingle(quantityFullname, c, units)
        } catch (error) {
          setError(`Invalid value for this metainfo. Please check your syntax.`)
          return
        }
        valid = true
      }
    }

    // Check if filter is locked
    if (filtersLocked[quantityFullname]) {
      setError(`Cannot change the filter as it is locked in the current search context.`)
      return
    }

    if (valid) {
      // Submit to search context on successful validation.
      setFilter([quantityFullname, old => {
        const multiple = filterData[quantityFullname].multiple
        return (isNil(old) || !multiple) ? queryValue : new Set([...old, ...queryValue])
      }])
      setInputValue('')
      setOpen(false)
    } else {
      setError(`Invalid query`)
    }
  }, [inputValue, filtersLocked, filters, units, setFilter, filterData])

  // Handle clear button
  const handleClose = useCallback(() => {
    setInputValue('')
    setOpen(false)
  }, [])

  const handleHighlight = useCallback((event, value, reason) => {
    setHighlighted(value)
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
    value = value?.trim()
    if (!value) {
      setSuggestionInput('')
      setOpen(false)
      return
    } else {
      setOpen(true)
    }
    if (reason !== 'input') {
      setSuggestionInput('')
      setOpen(false)
    }

    // If some input is given, and the quantity supports suggestions, we use
    // input suggester to suggest values. If the input is prefixed with a proper
    // quantity name and an equals-sign, we extract the quantity name and the
    // typed input
    const split = value.split('=', 2)
    let quantityList = [...quantitySetSuggestable]
    if (split.length === 2) {
      const quantityName = split[0].trim()
      const quantityFullname = filterFullnames[quantityName]
      if (quantitySetSuggestable.has(quantityName)) {
        quantityList = [quantityName]
        value = split[1].trim()
      } else if (quantitySetSuggestable.has(quantityFullname)) {
        quantityList = [quantityFullname]
        value = split[1].trim()
      }
    }
    setQuantityList(quantityList)
    setSuggestionInput(value)
  }, [quantitySetSuggestable])

  // This determines the order: notice that items should be sorted by group
  // first in order for the grouping to work correctly.
  const options = useMemo(() => {
    for (let suggestion of suggestions) {
      suggestion.value = `${suggestion.category}=${suggestion.value}`
    }
    return suggestions.concat(quantitySuggestions)
  }, [quantitySuggestions, suggestions])

  return <Paper className={clsx(className, styles.root)}>
    <Autocomplete
      className={styles.input}
      freeSolo
      clearOnBlur={false}
      inputValue={inputValue}
      value={null}
      open={open}
      onOpen={() => { if (inputValue.trim() !== '') { setOpen(true) } }}
      onClose={() => setOpen(false)}
      fullWidth
      disableClearable
      PaperComponent={CustomPaper}
      classes={{endAdornment: styles.endAdornment}}
      groupBy={(option) => option.category}
      filterOptions={filterOptions}
      options={options}
      onInputChange={handleInputChange}
      onHighlightChange={handleHighlight}
      getOptionLabel={option => option.value}
      getOptionSelected={(option, value) => false}
      renderInput={(params) => (
        <TextField
          {...params}
          className={styles.textField}
          variant="outlined"
          placeholder="Start typing a query or a keyword to get relevant suggestions."
          label={error || undefined}
          error={!!error}
          onKeyDown={handleEnter}
          InputLabelProps={{ shrink: true }}
          InputProps={{
            ...params.InputProps,
            classes: {
              notchedOutline: styles.notchedOutline
            },
            startAdornment: <Tooltip title="Add filter">
              <IconButton onClick={handleSubmit} className={styles.iconButton} aria-label="search">
                <SearchIcon />
              </IconButton>
            </Tooltip>,
            endAdornment: (<>
              {loading ? <CircularProgress color="inherit" size={20} /> : null}
              {(inputValue?.length || null) && <>
                <Tooltip title="Clear">
                  <IconButton onClick={handleClose} className={styles.iconButton} aria-label="clear">
                    <CloseIcon />
                  </IconButton>
                </Tooltip>
              </>}
            </>)
          }}
        />
      )}
    />
  </Paper>
})

SearchBar.propTypes = {
  className: PropTypes.string
}

export default SearchBar
