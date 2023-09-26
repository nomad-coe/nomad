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
import assert from 'assert'
import { isNil, has } from 'lodash'
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
import { useSearchContext } from './SearchContext'
import { searchQuantities } from '../../config'
import { quantityNameSearch } from './FilterRegistry'

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

const numericTypes = new Set([DType.Timestamp, DType.Int, DType.Float])

// Customized paper component for the autocompletion options
export const CustomPaper = (props) => {
  return <Paper elevation={3} {...props} />
}

export const useStyles = makeStyles(theme => ({
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
 * This component shows a search bar with autocomplete functionality.
 */
const SearchBar = React.memo(({
  quantities,
  className
}) => {
  const styles = useStyles()
  const units = useUnits()
  const {
    filters,
    filterData,
    filterFullnames,
    useSetFilters,
    useFiltersLocked,
    useParseQuery
  } = useSearchContext()
  const [inputValue, setInputValue] = useState('')
  const [suggestionInput, setSuggestionInput] = useState('')
  const [highlighted, setHighlighted] = useState({value: ''})
  const [open, setOpen] = useState(false)
  const [error, setError] = useState(false)
  const filtersLocked = useFiltersLocked()
  const setFilter = useSetFilters()
  const parseQuery = useParseQuery()

  const [quantitiesAll, quantitiesAllSet, quantitiesSuggestable] = useMemo(() => {
    let quantitySetSuggestable, quantityList
    // Custom list of quantities. They are validated against the current search
    // context.
    if (quantities) {
      for (const quantity of quantities) {
        assert(
          quantity.name === quantityNameSearch || has(filterData, quantity.name),
          `Quantity ${quantity.name} does not exist in the current search context.`
        )
      }
      quantityList = [...quantities]
      quantitySetSuggestable = new Set(
        quantityList
          .map(q => q.name)
          .filter(name => name === quantityNameSearch || searchQuantities[name]?.suggestion)
      )
    // Default quantities to use. Certain quantities are prioritized, others
    // come in alphabetical order.
    } else {
      quantitySetSuggestable = new Set([
        'results.material.elements',
        'results.material.chemical_formula_hill',
        'results.material.chemical_formula_anonymous',
        'results.material.structural_type',
        'results.material.symmetry.structure_name',
        'results.method.simulation.program_name',
        'authors.name'
      ])
      const filterList = [...filters]
      filterList
        .filter(q => !quantitySetSuggestable.has(q) && searchQuantities[q]?.suggestion)
        .forEach(q => quantitySetSuggestable.add(q))
      quantityList = filterList.map((name) => ({name, size: 5}))
      // The list of available quantity names is provided at the very
      // bottom.
      quantitySetSuggestable.add(quantityNameSearch)
      quantityList.push({name: quantityNameSearch})
    }
    const quantitySet = new Set(quantityList.map((q) => q.name))
    return [quantityList, quantitySet, quantitySetSuggestable]
  }, [quantities, filterData, filters])

  const [suggestionNames, setSuggestionNames] = useState(quantitiesSuggestable)
  const suggestionQuantities = useMemo(() => {
    return quantitiesAll.filter((q) => suggestionNames.has(q.name))
  }, [quantitiesAll, suggestionNames])
  const [suggestions, loading] = useSuggestions(suggestionQuantities, quantitiesAllSet, suggestionInput)

  // Used to check the validity of the given quantity name
  const checkMetainfo = useCallback((name) => {
    const fullName = filterFullnames[name] || name
    if (!quantitiesAllSet.has(fullName)) {
      return [undefined, `The quantity '${name}' is not supported`]
    }
    if (filterData[fullName].section) {
      return [fullName, `Cannot target metainfo sections`]
    }
    return [fullName, undefined]
  }, [quantitiesAllSet, filterData, filterFullnames])

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
      const [fullName, error] = checkMetainfo(quantityName)
      quantityFullname = fullName
      if (error) {
        setError(error)
        return
      }
      try {
        queryValue = parseQuery(quantityFullname, equals[2], units)
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
        const [aFullname, aError] = checkMetainfo(a)
        const [bFullname, bError] = checkMetainfo(b)
        if (aError && bError) {
          if (!aFullname && !bFullname) {
            setError(`Unknown quantity name`)
          } else {
            setError(aError || bError)
          }
          return
        }
        quantityFullname = aError ? bFullname : aFullname
        const value = aError ? a : b
        const dtype = getDatatype(quantityFullname)
        if (!numericTypes.has(dtype)) {
          setError(`Cannot perform range query for a non-numeric quantity`)
          return
        }
        let quantityValue
        try {
          quantityValue = parseQuery(quantityFullname, value, units)
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
        const [fullName, error] = checkMetainfo(b)
        if (error) {
          setError(error)
          return
        }
        quantityFullname = fullName
        const dtype = getDatatype(quantityFullname)
        if (!numericTypes.has(dtype)) {
          setError(`Cannot perform range query for a non-numeric quantity.`)
          return
        }
        queryValue = {}
        try {
          queryValue[opMapReverse[op1]] = parseQuery(quantityFullname, a, units)
          queryValue[opMap[op2]] = parseQuery(quantityFullname, c, units)
        } catch (error) {
          setError(`Invalid value for this metainfo. Please check your syntax.`)
          return
        }
        valid = true
      }
    }

    // Check if filter is locked
    if (!isNil(filtersLocked[quantityFullname]) && filterData[quantityFullname].global) {
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
  }, [inputValue, checkMetainfo, units, setFilter, filterData, parseQuery, filtersLocked])

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
      if (open && highlighted?.text) {
        setInputValue(highlighted.text)
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
    let quantitySet = quantitiesSuggestable
    if (split.length === 2) {
      const quantityName = split[0].trim()
      const quantityFullname = filterFullnames[quantityName]
      if (quantitiesSuggestable.has(quantityName)) {
        quantitySet = new Set([quantityName])
        value = split[1].trim()
      } else if (quantitiesSuggestable.has(quantityFullname)) {
        quantitySet = new Set([quantityFullname])
        value = split[1].trim()
      }
    }
    setSuggestionNames(quantitySet)
    setSuggestionInput(value)
  }, [quantitiesSuggestable, filterFullnames])

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
      options={suggestions}
      onInputChange={handleInputChange}
      onHighlightChange={handleHighlight}
      getOptionLabel={option => option.text}
      // Notice that we need to override the default filterOptions as it is
      // performing unwanted filtering.
      filterOptions={(options) => options}
      renderInput={(params) => (
        <TextField
          {...params}
          className={styles.textField}
          variant="outlined"
          placeholder="Type your query or keyword here"
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
 // List of quantities that the search targets. The suggestions are returned in
 // the same order.
  quantities: PropTypes.arrayOf(PropTypes.shape({
    name: PropTypes.string,
    size: PropTypes.number
  })),
  className: PropTypes.string
}

export default SearchBar
