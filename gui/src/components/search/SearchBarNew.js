import React, { useRef, useState, useContext, useCallback, useMemo } from 'react'
import {searchContext} from './SearchContext'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TextField from '@material-ui/core/TextField'
import { CircularProgress } from '@material-ui/core'
import * as searchQuantities from '../../searchQuantities.json'
import { apiContext } from '../api'

export default function SearchBar() {
  const suggestionsTimerRef = useRef(null)
  const {response: {statistics, pagination}, domain, query, setQuery} = useContext(searchContext)
  const defaultOptions = useMemo(() => {
    return Object.keys(searchQuantities)
      .map(quantity => searchQuantities[quantity].name)
      .filter(quantity => !quantity.includes('.') || quantity.startsWith(domain.key + '.'))
  }, [domain.key])

  const [open, setOpen] = useState(false)
  const [options, setOptions] = useState(defaultOptions)
  const [loading, setLoading] = useState(false)
  const [inputValue, setInputValue] = useState('')

  const {api} = useContext(apiContext)

  const autocompleteValue = Object.keys(query).map(quantity => `${quantity}=${query[quantity]}`)

  let helperText = ''
  if (pagination && statistics) {
    if (pagination.total === 0) {
      helperText = <span>There are no more entries matching your criteria.</span>
    } else {
      helperText = <span>
        There {pagination.total === 1 ? 'is' : 'are'} {Object.keys(domain.searchMetrics).filter(key => statistics.total.all[key]).map(key => {
          return <span key={key}>
            {domain.searchMetrics[key].renderResultString(statistics.total.all[key])}
          </span>
        })}{Object.keys(query).length ? ' left' : ''}.
      </span>
    }
  }

  const filterOptions = useCallback((options, params) => {
    const [quantity, value] = params.inputValue.split('=')
    const filteredOptions = options.filter(option => {
      const [optionQuantity, optionValue] = option.split('=')
      if (!value) {
        return optionQuantity.includes(quantity) || optionQuantity === quantity
      } else {
        return optionValue.includes(value) || optionValue === value
      }
    })
    return filteredOptions
  }, [])

  const loadOptions = useCallback((quantity, value) => {
    if (suggestionsTimerRef.current !== null) {
      clearTimeout(suggestionsTimerRef.current)
    }
    suggestionsTimerRef.current = setTimeout(() => {
      setLoading(true)
      api.suggestions_search(quantity, query, value, 20, true)
        .then(response => {
          setLoading(false)
          const options = response.suggestions.map(value => `${quantity}=${value}`)
          setOptions(options)
          setOpen(true)
        })
        .catch(() => {
          setLoading(false)
        })
    }, 200)
  }, [api, suggestionsTimerRef])

  const handleInputChange = useCallback((event, value, reason) => {
    if (reason === 'input') {
      setInputValue(value)
      const [quantity, quantityValue] = value.split('=')

      if (searchQuantities[quantity]) {
        loadOptions(quantity, quantityValue)
      } else {
        setOptions(defaultOptions)
      }
    }
  }, [loadOptions])

  const handleChange = (event, entries) => {
    const newQuery = entries.reduce((query, entry) => {
      if (entry) {
        const [quantity, value] = entry.split('=')
        if (query[quantity]) {
          if (searchQuantities[quantity].many) {
            if (Array.isArray(query[quantity])) {
              query[quantity].push(value)
            } else {
              query[quantity] = [query[quantity], value]
            }
          } else {
            query[quantity] = value
          }
        } else {
          query[quantity] = value
        }
      }
      return query
    }, {})
    setQuery(newQuery, true)

    if (entries.length !== 0) {
      const entry = entries[entries.length - 1]
      const [quantity, value] = entry.split('=')
      if (value) {
        setInputValue('')
      } else {
        setInputValue(`${entry}=`)
        loadOptions(quantity)
      }
    }
  }

  React.useEffect(() => {
    if (!open) {
      setOptions(defaultOptions)
    }
  }, [open])

  return <Autocomplete
    multiple
    freeSolo
    inputValue={inputValue}
    value={autocompleteValue}
    limitTags={4}
    id='search-bar'
    open={open}
    onOpen={() => {
      setOpen(true)
    }}
    onClose={() => {
      setOpen(false)
    }}
    onChange={handleChange}
    onInputChange={handleInputChange}
    getOptionSelected={(option, value) => option === value}
    options={options}
    loading={loading}
    filterOptions={filterOptions}
    renderInput={(params) => (
      <TextField
        {...params}
        helperText={helperText}
        label='Search with quantity=value'
        variant='outlined'
        InputProps={{
          ...params.InputProps,
          endAdornment: (
            <React.Fragment>
              {loading ? <CircularProgress color='inherit' size={20} /> : null}
              {params.InputProps.endAdornment}
            </React.Fragment>
          )
        }}
      />
    )}
  />
}
