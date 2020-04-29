import React from 'react'
import {searchContext} from './SearchContext'
import Autocomplete from '@material-ui/lab/Autocomplete'
import TextField from '@material-ui/core/TextField'
import { CircularProgress } from '@material-ui/core'
import * as searchQuantities from '../../searchQuantities.json'
import { apiContext } from '../api'

const defaultOptions = Object.keys(searchQuantities).map(quantity => searchQuantities[quantity].name)

export default function SearchBar() {
  const [open, setOpen] = React.useState(false)
  const [options, setOptions] = React.useState(defaultOptions)
  const [loading, setLoading] = React.useState(false)
  const [inputValue, setInputValue] = React.useState('')
  const {response: {statistics, pagination}, domain, query, setQuery} = React.useContext(searchContext)
  const {api} = React.useContext(apiContext)

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

  const loadValues = (quantity, value) => {
    setLoading(true)
    const size = searchQuantities[quantity].statistic_size || 20
    api.quantity_search(quantity, query, size, true)
      .then(response => {
        setLoading(false)
        const options = Object.keys(response.quantity.values).map(value => `${quantity}=${value}`)
        setOptions(options)
        setOpen(true)
      })
      .catch(() => {
        setLoading(false)
      })
  }

  const handleInputChange = (event, value, reason) => {
    if (reason === 'input') {
      setInputValue(value)
      const [quantity, quantityValue] = value.split('=')
      if (!quantityValue) {
        if (searchQuantities[quantity]) {
          loadValues(quantity, quantityValue)
        } else {
          setOptions(defaultOptions)
        }
      }
    }
  }

  const handleChange = (event, entries, reason) => {
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
        loadValues(quantity)
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
    getOptionLabel={(option) => option}
    options={options}
    loading={loading}
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
