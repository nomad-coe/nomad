import React, { useRef, useState, useContext, useCallback, useMemo } from 'react'
import {searchContext} from './SearchContext'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TextField from '@material-ui/core/TextField'
import { CircularProgress, InputAdornment, Button, Tooltip } from '@material-ui/core'
import searchQuantities from '../../searchQuantities'
import { apiContext } from '../api'

/**
 * A few helper functions related to format and analyse suggested options
 */
const Options = {
  split: (suggestion) => suggestion.split('='),
  join: (quantity, value) => `${quantity}=${value}`,
  splitForCompare: (suggestion) => {
    const [quantity, value] = suggestion.split('=')
    return [quantity ? quantity.toLowerCase() : '', value ? value.toLowerCase() : '']
  }
}

/**
 * This searchbar component shows a searchbar with autocomplete functionality. The
 * searchbar also includes a status line about the current results. It uses the
 * search context to manipulate the current query and display results. It does its on
 * API calls to provide autocomplete suggestion options.
 */
export default function SearchBar() {
  const suggestionsTimerRef = useRef(null)
  const {response: {statistics, pagination, error}, domain, query, apiQuery, setQuery} = useContext(searchContext)
  const defaultOptions = useMemo(() => {
    return Object.keys(searchQuantities)
      .map(quantity => searchQuantities[quantity].name)
      .filter(quantity => !quantity.includes('.') || quantity.startsWith(domain.key + '.'))
  }, [domain.key])

  const [open, setOpen] = useState(false)
  const [options, setOptions] = useState(defaultOptions)
  const [loading, setLoading] = useState(false)
  const [inputValue, setInputValue] = useState('')
  const [searchType, setSearchType] = useState('nomad')

  const {api} = useContext(apiContext)

  const autocompleteValue = Object.keys(query).map(quantity => Options.join(quantity, query[quantity]))

  const handleSearchTypeClicked = useCallback(() => {
    if (searchType === 'nomad') {
      setSearchType('optimade')
      // handleChange(null, [])
    } else {
      setSearchType('nomad')
      // handleChange(null, [])
    }
  }, [searchType, setSearchType])

  const handleOptimadeEntered = useCallback(query => {
    setQuery({'dft.optimade': query})
  })

  let helperText = ''
  if (error) {
    helperText = '' + (error.apiMessage || error)
  } else if (pagination && statistics) {
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
    let [quantity, value] = Options.splitForCompare(params.inputValue)
    const filteredOptions = options.filter(option => {
      let [optionQuantity, optionValue] = Options.splitForCompare(option)
      if (!value) {
        return optionQuantity && (optionQuantity.includes(quantity) || optionQuantity === quantity)
      } else {
        return optionValue.includes(value) || optionValue === value
      }
    })
    return filteredOptions
  }, [])

  const loadOptions = useCallback((quantity, value) => {
    const size = searchQuantities[quantity].statistic_size

    if (suggestionsTimerRef.current !== null) {
      clearTimeout(suggestionsTimerRef.current)
    }
    if (loading) {
      return
    }
    suggestionsTimerRef.current = setTimeout(() => {
      setLoading(true)
      api.suggestions_search(quantity, apiQuery, size ? null : value, size || 20, true)
        .then(response => {
          setLoading(false)
          const options = response.suggestions.map(value => Options.join(quantity, value))
          setOptions(options)
          setOpen(true)
        })
        .catch(() => {
          setLoading(false)
        })
    }, 200)
  }, [api, suggestionsTimerRef, apiQuery])

  const handleInputChange = useCallback((event, value, reason) => {
    if (reason === 'input') {
      setInputValue(value)
      const [quantity, quantityValue] = Options.split(value)

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
        const [quantity, value] = Options.split(entry)
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
      const [quantity, value] = Options.split(entry)
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

  const commonTextFieldProps = params => ({
    error: !!error,
    helperText: helperText,
    variant: 'outlined',
    fullWidth: true,
    ...params
  })

  const commonInputProps = (params) => ({
    ...params,
    startAdornment: (
      <React.Fragment>
        <InputAdornment position="start">
          <Tooltip title="Switch between NOMAD's quantity=value search and the Optimade filter language.">
            <Button onClick={handleSearchTypeClicked}size="small">{searchType}</Button>
          </Tooltip>
        </InputAdornment>
        {params.startAdornment}
      </React.Fragment>
    )
  })

  if (searchType === 'nomad') {
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
          {...commonTextFieldProps(params)}
          label={searchType === 'nomad' ? 'Search with quantity=value' : 'Search with Optimade filter language'}
          InputProps={{
            ...commonInputProps(params.InputProps),
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
  } else {
    return <TextField
      {...commonTextFieldProps({})}
      label={searchType === 'nomad' ? 'Search with quantity=value' : 'Search with Optimade filter language'}
      InputProps={{
        ...commonInputProps({})
      }}
      defaultValue={query['dft.optimade'] || ''}
      onKeyPress={(ev) => {
        console.log(`Pressed keyCode ${ev.key}`)
        if (ev.key === 'Enter') {
          handleOptimadeEntered(ev.target.value)
          ev.preventDefault()
        }
      }}
    />
  }
}
