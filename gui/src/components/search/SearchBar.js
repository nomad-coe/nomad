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
import React, { useRef, useState, useContext, useCallback, useMemo } from 'react'
import {searchContext} from './SearchContext'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TextField from '@material-ui/core/TextField'
import { CircularProgress, InputAdornment, Button, Tooltip } from '@material-ui/core'
import searchQuantities from '../../searchQuantities'
import { apiContext } from '../api'
import { defsByName as metainfoDefs } from '../archive/metainfo'
import { domains } from '../domains'

const metainfoOptions = []

const quantitiesWithAlternativeOptions = {
  calc_id: () => [],
  upload_id: () => [],
  calc_hash: () => [],
  'dft.quantities': () => {
    if (metainfoOptions.length === 0) {
      metainfoOptions.push(...Object.keys(metainfoDefs)
        .filter(name => !name.startsWith('x_'))
        .map(name => ({
          domain: 'dft',
          quantity: 'dft.quantities',
          value: name
        })))
    }
    return metainfoOptions
  }
}

// We need to treat dft. and encyclopedia. special. Usually all dft domain pieces
// are prefixed dft., but the encycloepdia is top-level and also a dft. specific
// quantity. These to functions remove and add the dft./encyclopedia. prefixes accordingly.
function getDomainOfQuantity(quantity) {
  if (!quantity.includes('.')) {
    return null
  }
  const firstSegment = quantity.split('.')[0]
  if (firstSegment === 'encyclopedia') {
    return 'dft'
  }
  return firstSegment
}

function addDomainToQuantity(shortenedQuantityName, domainKey) {
  if (!searchQuantities[shortenedQuantityName]) {
    shortenedQuantityName = domainKey + '.' + shortenedQuantityName
    if (!searchQuantities[shortenedQuantityName]) {
      shortenedQuantityName = 'encyclopedia.' + shortenedQuantityName.slice(4)
    }
  }
  return shortenedQuantityName
}

/**
 * This searchbar component shows a searchbar with autocomplete functionality. The
 * searchbar also includes a status line about the current results. It uses the
 * search context to manipulate the current query and display results. It does its on
 * API calls to provide autocomplete suggestion options.
 */
export default function SearchBar() {
  const currentLoadOptionsConfigRef = useRef({
    timer: null,
    latestOption: null,
    requestedOption: null
  })
  const {response: {statistics, pagination, error}, domain, query, apiQuery, setQuery} = useContext(searchContext)
  const defaultOptions = useMemo(() => {
    return Object.keys(searchQuantities)
      .map(quantity => ({
        quantity: quantity,
        domain: getDomainOfQuantity(quantity)
      }))
      .filter(option => !option.domain || option.domain === domain.key)
  }, [domain.key])

  const [open, setOpen] = useState(false)
  const [options, setOptions] = useState(defaultOptions)
  const [loading, setLoading] = useState(false)
  const [inputValue, setInputValue] = useState('')
  const [searchType, setSearchType] = useState('nomad')

  const {api} = useContext(apiContext)

  const autocompleteValue = Object.keys(query).map(quantity => ({
    quantity: quantity,
    domain: quantity.includes('.') ? quantity.split('.')[0] : null,
    value: query[quantity]
  }))

  const handleSearchTypeClicked = useCallback(() => {
    if (searchType === 'nomad') {
      setSearchType('optimade')
    } else {
      setSearchType('nomad')
    }
  }, [searchType, setSearchType])

  const handleOptimadeEntered = useCallback(query => {
    setQuery({'dft.optimade': query})
  }, [setQuery])

  let helperText = ''
  if (error) {
    helperText = '' + (error.apiMessage || error)
  } else if (pagination && statistics) {
    if (pagination.total === 0) {
      helperText = <span>There are no more entries matching your criteria.</span>
    } else {
      helperText = <span>
        There {pagination.total === 1 ? 'is' : 'are'} {
          Object.keys(domain.searchMetrics).filter(key => statistics.total.all[key]).map(key => {
            return <span key={key}>
              {domain.searchMetrics[key].renderResultString(statistics.total.all[key])}
            </span>
          })
        }{Object.keys(query).length ? ' left' : ''}.
      </span>
    }
  }

  const loadOptions = useCallback(option => {
    const config = currentLoadOptionsConfigRef.current
    config.latestOption = option

    if (config.timer !== null) {
      clearTimeout(config.timer)
    }
    if (loading) {
      return
    }
    config.timer = setTimeout(() => {
      config.requestedOption = option

      const alternativeOptions = quantitiesWithAlternativeOptions[option.quantity]
      if (alternativeOptions) {
        setOptions(alternativeOptions())
        return
      }

      const size = searchQuantities[option.quantity].statistic_size
      setLoading(true)
      api.suggestions_search(option.quantity, apiQuery, size ? null : option.value, size || 20, true)
        .then(response => {
          setLoading(false)
          if (!config.latestOption || config.requestedOption.quantity !== config.latestOption.quantity) {
            // don't do anything if quantity has changed in the meantime
            return
          }
          const options = response.suggestions.map(value => ({
            quantity: option.quantity,
            domain: option.domain,
            value: value
          }))
          setOptions(options)
          setOpen(true)
        })
        .catch(() => {
          setLoading(false)
        })
    }, 200)
  }, [api, currentLoadOptionsConfigRef, apiQuery, loading, setLoading])

  const getOptionLabel = useCallback(option => {
    if (option.quantity === 'from_time' || option.quantity === 'until_time') {
      if (option.value) {
        return `${option.quantity.replace('_time', '')}=${option.value.substring(0, 10)}`
      }
    }

    let label = option.quantity + '='
    if (option.value) {
      if (Array.isArray(option.value)) {
        label += option.value.join(',')
      } else {
        label += option.value
      }
    }
    return label.substring(label.indexOf('.') + 1)
  }, [])

  const parseOption = useCallback(input => {
    const [inputQuantity, inputValue] = input.split('=')

    const quantity = addDomainToQuantity(inputQuantity, domain.key)
    let value = inputValue
    if (value && searchQuantities[quantity] && searchQuantities[quantity].many) {
      value = value.split(',').map(item => item.trim())
    }
    return {
      inputQuantity: inputQuantity,
      inputValue: inputValue,
      domain: inputQuantity.includes('.') ? inputQuantity.split('.')[0] : null,
      quantity: searchQuantities[quantity] ? quantity : null,
      value: value
    }
  }, [domain.key])

  const filterOptions = useCallback((options, params) => {
    const inputOption = parseOption(params.inputValue)
    const filteredOptions = options.filter(option => {
      if (!inputOption.quantity) {
        return option.quantity.includes(
          inputOption.inputQuantity) && (option.domain === domain.key || !option.domain)
      }
      if (option.quantity !== inputOption.quantity) {
        return false
      }

      if (!inputOption.value) {
        return true
      }

      const matches = option.value &&
        inputOption.inputValue &&
        option.value.toLowerCase().includes(inputOption.inputValue.toLowerCase())
      if (matches) {
        if (option.value === inputOption.inputValue) {
          inputOption.exists |= true
        }
        return true
      }

      return false
    })

    // Add the value as option, even if it does not exist to allow search for missing,
    // faulty, or not yet loaded options
    if (inputOption.quantity && !inputOption.exists) {
      filteredOptions.push(inputOption)
    }

    return filteredOptions
  }, [domain.key, parseOption])

  const handleInputChange = useCallback((event, value, reason) => {
    if (reason === 'input') {
      setInputValue(value)
      const inputOption = parseOption(value)
      if (inputOption.quantity) {
        loadOptions(inputOption)
      } else {
        setOptions(defaultOptions)
      }
    }
  }, [loadOptions, defaultOptions, parseOption])

  const handleChange = (event, entries) => {
    currentLoadOptionsConfigRef.current.latestOption = null

    entries = entries.map(entry => {
      if (typeof entry === 'string') {
        return parseOption(entry)
      } else {
        return entry
      }
    })

    const newQuery = entries.reduce((query, entry) => {
      if (entry) {
        if (query[entry.quantity]) {
          if (searchQuantities[entry.quantity].many) {
            if (Array.isArray(query[entry.quantity])) {
              query[entry.quantity].push(entry.value)
            } else {
              query[entry.quantity] = [query[entry.quantity], entry.value]
            }
          } else {
            query[entry.quantity] = entry.value
          }
        } else {
          query[entry.quantity] = entry.value
        }
      }
      return query
    }, {})
    setQuery(newQuery, true)

    if (entries.length !== 0) {
      const entry = entries[entries.length - 1]
      if (entry.value) {
        setInputValue('')
      } else {
        setInputValue(getOptionLabel(entry))
        loadOptions(entry)
      }
    }
  }

  React.useEffect(() => {
    if (!open) {
      setOptions(defaultOptions)
    }
  }, [open, defaultOptions])

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
        {domain === domains.dft &&
        <InputAdornment position="start">
          <Tooltip title="Switch between NOMAD's quantity=value search and the Optimade filter language.">
            <Button onClick={handleSearchTypeClicked}size="small">{searchType}</Button>
          </Tooltip>
        </InputAdornment>}
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
      getOptionSelected={(option, inputOption) => {
        return inputOption.quantity === option.quantity && inputOption.value === option.value
      }}
      getOptionLabel={getOptionLabel}
      options={options}
      loading={loading}
      filterOptions={filterOptions}
      // handleHomeEndKeys
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
        if (ev.key === 'Enter') {
          handleOptimadeEntered(ev.target.value)
          ev.preventDefault()
        }
      }}
    />
  }
}
