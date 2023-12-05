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
import { makeStyles } from '@material-ui/core/styles'
import SearchIcon from '@material-ui/icons/Search'
import { Paper, Tooltip } from '@material-ui/core'
import IconButton from '@material-ui/core/IconButton'
import { DType } from '../../utils'
import { useSuggestions } from '../../hooks'
import { useSearchContext } from './SearchContext'
import { quantityNameSearch } from './FilterRegistry'
import { MetainfoOption, ListboxMetainfo, renderGroup } from './input/InputMetainfo'
import { InputText } from './input/InputText'

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
  const {
    filters,
    filterData,
    filterFullnames,
    useUpdateFilter,
    useFiltersLocked,
    useParseQuery
  } = useSearchContext()
  const [inputValue, setInputValue] = useState('')
  const [suggestionInput, setSuggestionInput] = useState('')
  const [error, setError] = useState(false)
  const filtersLocked = useFiltersLocked()
  const updateFilter = useUpdateFilter()
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
          .filter(name => name === quantityNameSearch || filterData[name]?.suggestion)
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
      filterList.forEach(q => quantitySetSuggestable.add(q))
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
  const [suggestions, loading] = useSuggestions(suggestionQuantities, quantitiesAllSet, suggestionInput, filterData)
  const {options, keys} = useMemo(() => {
    const keys = []
    const options = {}
    for (const suggestion of suggestions) {
      let filter, key, primary
      if (suggestion.category === quantityNameSearch) {
        filter = filterData[suggestion.value]
        key = suggestion.value
        primary = filter?.quantity || suggestion.value
      } else {
        key = `${suggestion.category}=${suggestion.value}`
        primary = key
      }
      keys.push(key)
      options[key] = {
        key: key,
        group: suggestion.category,
        dtype: filter?.dtype,
        schema: filter?.schema,
        primary: primary,
        description: filter?.description
      }
    }
    return {options, keys}
  }, [suggestions, filterData])

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
  const handleSubmit = useCallback((value) => {
    if (!value?.trim()?.length) {
      return
    }
    const reString = '[^\\s=<>](?:[^=<>]*[^\\s=<>])?'
    const op = '(?:<|>)=?'
    let valid = false
    let quantityFullname
    let queryValue
    let comparison = true

    // Presence query
    const presence = value.match(new RegExp(`^\\s*(${reString})\\s*=\\s*\\*\\s*$`))
    if (presence) {
      quantityFullname = `quantities`
      queryValue = parseQuery(quantityFullname, presence[1])
      valid = true
    }

    // Equality query
    if (!valid) {
      const equals = value.match(new RegExp(`^\\s*(${reString})\\s*=\\s*(${reString})\\s*$`))
      if (equals) {
        const quantityName = equals[1]
        const [fullName, error] = checkMetainfo(quantityName)
        quantityFullname = fullName
        if (error) {
          setError(error)
          return
        }
        try {
          queryValue = parseQuery(quantityFullname, equals[2])
        } catch (error) {
          setError(`Invalid value for this metainfo. Please check your syntax.`)
          return
        }
        comparison = false
        valid = true
      }
    }

    // Simple LTE/GTE query
    if (!valid) {
      const ltegte = value.match(new RegExp(`^\\s*(${reString})\\s*(${op})\\s*(${reString})\\s*$`))
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
        const dtype = filterData[quantityFullname].dtype
        if (!numericTypes.has(dtype)) {
          setError(`Cannot perform range query for a non-numeric quantity`)
          return
        }
        let quantityValue
        try {
          quantityValue = parseQuery(quantityFullname, value, undefined, false)
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
      const ltegteSandwich = value.match(new RegExp(`^\\s*(${reString})\\s*(${op})\\s*(${reString})\\s*(${op})\\s*(${reString})\\s*$`))
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
        const dtype = filterData[quantityFullname].dtype
        if (!numericTypes.has(dtype)) {
          setError(`Cannot perform range query for a non-numeric quantity.`)
          return
        }
        queryValue = {}
        try {
          queryValue[opMapReverse[op1]] = parseQuery(quantityFullname, a, undefined, false)
          queryValue[opMap[op2]] = parseQuery(quantityFullname, c, undefined, false)
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

    // Submit to search context on successful validation.
    if (valid) {
      updateFilter([quantityFullname, old => {
        const multiple = filterData[quantityFullname].multiple
        return (comparison || isNil(old) || !multiple)
          ? queryValue
          : new Set([...old, ...queryValue])
      }])
      setInputValue('')
    } else {
      setError(`Invalid query`)
    }
  }, [checkMetainfo, updateFilter, filterData, parseQuery, filtersLocked])

  const handleSelect = useCallback((key) => {
    setInputValue(key)
  }, [])

  // Handle typing events. After a debounce time has expired, a list of
  // suggestion will be retrieved if they are available for this metainfo and
  // the input is deemed meaningful.
  const handleInputChange = useCallback((value) => {
    setError(error => error ? undefined : null)
    setInputValue(value)
    value = value?.trim()
    if (!value) {
      setSuggestionInput('')
      return
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
    <InputText
      value={inputValue || null}
      onChange={handleInputChange}
      onSelect={handleSelect}
      onAccept={handleSubmit}
      suggestions={keys}
      ListboxComponent={ListboxMetainfo}
      TextFieldProps={{
        variant: 'outlined',
        placeholder: 'Type your query or keyword here',
        label: error || undefined,
        error: !!error,
        InputLabelProps: { shrink: true },
        size: "medium"
      }}
      InputProps={{
        classes: {
          notchedOutline: styles.notchedOutline
        },
        startAdornment: <Tooltip title="Add filter">
          <IconButton onClick={() => handleSubmit(inputValue)} className={styles.iconButton} aria-label="search">
            <SearchIcon />
          </IconButton>
        </Tooltip>
      }}
      groupBy={(key) => options?.[key]?.group}
      renderGroup={renderGroup}
      getOptionLabel={option => option}
      filterOptions={(options) => options}
      loading={loading}
      renderOption={(id) => <MetainfoOption id={id} options={options} />}
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
