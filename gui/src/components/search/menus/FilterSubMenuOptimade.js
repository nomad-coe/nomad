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
import React, {useCallback, useContext, useEffect, useMemo, useState} from 'react'
import PropTypes from 'prop-types'
import {FilterSubMenu, filterMenuContext} from './FilterMenu'
import {InputCheckboxValue} from '../input/InputCheckbox'
import {Box, TextField} from "@material-ui/core"
import {useSearchContext} from "../SearchContext"
import AutoComplete from "@material-ui/lab/Autocomplete"
import {useApi} from '../../api'
import {debounce} from 'lodash'
import {useErrors} from '../../errors'

const FilterSubMenuOptimade = React.memo(({
  value,
  ...rest
}) => {
  const [suggestions, setSuggestions] = useState([['']])
  const [identifiers, setIdentifiers] = useState([])
  const [values, setValues] = useState([{value: '', valid: true, msg: ''}])
  const {selected, open} = useContext(filterMenuContext)
  const {useFilterState} = useSearchContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const visible = open && value === selected
  const [optimadeFilters, setOptimadeFilters] = useFilterState("optimade_filter")

  useEffect(() => {
    const requestBody = {
      exclude: ['atoms', 'only_atoms', 'files', 'quantities', 'dft.quantities', 'dft.labels', 'dft.geometries'],
      owner: 'public'
    }
    api.post(`entries/query`, requestBody).then(response => {
      const optimade = response?.data?.[0]?.optimade
      const keys = optimade && Object.keys(optimade)
      setIdentifiers(keys || [])
    }).catch(error => {
      raiseError(error)
      setIdentifiers([])
    })
  }, [api, raiseError])

  const renderSuggestion = useCallback((command, options, suggestions) => {
    options.filter((item, pos) => options.indexOf(item) === pos).forEach(option => {
      if (option === 'IDENTIFIER') {
        identifiers.forEach(identifier => {
          suggestions.push([command, identifier].join(' '))
        })
      } else if (['OPERATOR', 'ESCAPED_STRING', 'SIGNED_FLOAT', 'SIGNED_INT'].includes(option)) {
        suggestions.push(command.concat(' <', option, '>'))
      } else {
        suggestions.push([command, option].join(' '))
      }
    })
  }, [identifiers])

  const getSuggestions = useCallback(async (filter, index) => {
    const requestBody = {
      exclude: ['atoms', 'only_atoms', 'files', 'quantities', 'dft.quantities', 'optimade', 'dft.labels', 'dft.geometries'],
      verify_only: true,
      owner: 'public',
      query: {
        optimade_filter: filter
      }
    }
    try {
      await api.post(`entries/query`, requestBody)
      setValues(oldValues => {
        const newValues = [...oldValues]
        newValues[index] = {value: oldValues[index].value, valid: true, msg: ''}
        return newValues
      })
      return {valid: true, suggestions: []}
    } catch (error) {
      if (error?.apiMessage) {
        const suggestions = []
        error?.apiMessage.forEach((err) => {
          const errorLocations = err?.loc
          errorLocations.forEach(location => {
            if (location?.[0] === 'optimade_filter') {
              let errorIndex = Number(err.msg.match(/col \d+/i)?.[0].match(/\d+/)?.[0]) - 1
              if (err.msg.match(/Unexpected end-of-input/i)) {
                errorIndex = filter.length
              }
              const commands = identifiers.filter(identifier => identifier.startsWith(filter))
              if (commands.length > 0) {
                commands.forEach(command => suggestions.push(command))
              }
              const options = Array.from(err.msg.matchAll(/\t\* \w+\n/g), m => m[0].slice(3, -1))
              const command = filter.slice(0, errorIndex).trim()
              renderSuggestion(command, options, suggestions)
              if (err.msg.match(/Semantic error/i)) {
                const msg = err.msg.match(/(?<=:\n\n).*/i)?.[0]
                setValues(oldValues => {
                  const newValues = [...oldValues]
                  newValues[index] = {value: oldValues[index].value, valid: oldValues[index].valid, msg: msg}
                  return newValues
                })
              }
            }
          })
        })
        return {valid: false, suggestions: suggestions}
      }
    }
  }, [api, identifiers, renderSuggestion])

  useEffect(() => {
    const newOptimadeFilters = optimadeFilters ? [...optimadeFilters, ''] : ['']
    const newValues = newOptimadeFilters.map(value => ({value: value, valid: !!value, msg: ''}))
    setValues(newValues)
    const prepareOptions = async () => {
      const allSuggestions = await Promise.all(newOptimadeFilters.map(async (optimadeFilter, index) => {
        const command = optimadeFilter && optimadeFilter.trim()
        const validation = await getSuggestions(command || '.', index)
        return [optimadeFilter].concat(validation.suggestions)
      }))
      setSuggestions(allSuggestions)
    }
    prepareOptions()
  }, [optimadeFilters, getSuggestions])

  const validate = useCallback(async (value, index) => {
    const validation = await getSuggestions((value && value.trim()) || '.', index)
    setSuggestions(suggestions => {
      const newSuggestions = [...suggestions]
      newSuggestions[index] = [value].concat(validation.suggestions)
      return newSuggestions
    })
    return validation
  }, [getSuggestions])

  const debouncedPrepareSuggestions = useMemo(() => (
    debounce(validate, 500)
  ), [validate])

  const setFilter = useCallback((value, index) => {
    setValues(oldValue => {
      const newValues = [...oldValue]
      if (!value) {
        index !== newValues.length - 1 && newValues.splice(index, 1)
      } else {
        newValues[index] = {value: value, valid: false, msg: ''}
      }
      return newValues
    })
    const newSuggestions = [...suggestions]
    if (!value) {
      index !== newSuggestions.length - 1 && newSuggestions.splice(index, 1)
    } else {
      newSuggestions[index] = [value]
    }
    setSuggestions(newSuggestions)
  }, [suggestions])

  const handleInputChange = useCallback((value, index) => {
    setFilter(value, index)
    debouncedPrepareSuggestions(value, index)
  }, [debouncedPrepareSuggestions, setFilter])

  const setOptimadeFilter = useCallback((index, value, valid = false) => {
    const newFilters = [...values].filter(value => !!value.value)
    if (newFilters.length === 0 || (!newFilters.every(value => !!value.valid) && !valid)) return
    const newOptimadeFilters = [...newFilters].map(value => value.value)
    if (!value) {
      newOptimadeFilters.splice(index, 1)
    } else {
      newOptimadeFilters[index] = value
    }
    setOptimadeFilters(new Set(newOptimadeFilters))
  }, [setOptimadeFilters, values])

  const handleChange = useCallback(async (event, value, index) => {
    if (values[index].value !== value) {
      setFilter(value, index)
      if (!values[index].value) return
    }
    if (!values[index].valid) {
      const validation = await validate(value, index)
      if (!validation.valid) return
      setOptimadeFilter(index, value, true)
    } else {
      setOptimadeFilter(index, value)
    }
  }, [setFilter, setOptimadeFilter, validate, values])

  return <FilterSubMenu
    value={value}
    actions={<InputCheckboxValue
      quantity="quantities"
      value="data"
      description="Search by optimade filters"
    />}
    {...rest}
  >
    {visible && values && values.length === suggestions.length && values.map((value, index) => {
      return <Box key={index} paddingLeft={2} paddingRight={2}>
        <AutoComplete
          options={suggestions[index]}
          style={{width: '100%'}}
          onChange={(event, value) => handleChange(event, value, index)}
          onKeyDown={(event) => (event.key === 'Enter' && handleChange(event, value.value, index))}
          value={value.value}
          getOptionSelected={(option, value) => option === value}
          getOptionLabel={option => (option && String(option)) || ''}
          renderInput={params => (
            <TextField
              {...params}
              variant='filled' label='Filter'
              placeholder='Type optimade filter' margin='normal' fullWidth size='small'
              onChange={event => handleInputChange(event.target.value, index)}
              error={!!value.msg}
              helperText={!!value.msg && value.msg}
            />
          )}
        />
      </Box>
    })}
  </FilterSubMenu>
})
FilterSubMenuOptimade.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuOptimade