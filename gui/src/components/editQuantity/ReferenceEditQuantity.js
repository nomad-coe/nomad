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
import React, { useCallback, useEffect, useMemo, useState, useRef } from 'react'
import PropTypes from 'prop-types'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { debounce } from 'lodash'
import { Autocomplete } from '@material-ui/lab'
import { TextField } from '@material-ui/core'
import { useEntryPageContext } from '../entry/EntryPageContext'
import { resolveRefAsync } from '../archive/metainfo'
import { ItemButton, useLane } from '../archive/Browser'
import { useBrowserAdaptorContext } from '../archive/ArchiveBrowser'
import { getFieldProps } from './StringEditQuantity'
import { isWaitingForUpdateTestId } from '../../utils'

const ReferenceEditQuantity = React.memo(function ReferenceEditQuantity(props) {
  const {archive} = useEntryPageContext()
  const {quantityDef, value, onChange, index} = props
  const [entry, setEntry] = useState(null)
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [inputValue, setInputValue] = useState('')
  const [suggestions, setSuggestions] = useState([])
  const [error, setError] = useState()
  const fetchedSuggestionsFor = useRef()
  const lane = useLane()
  const browserAdaptorContext = useBrowserAdaptorContext(archive)
  const context = lane?.adaptor?.context || browserAdaptorContext
  const referencedSectionQualifiedName = useMemo(() => {
    return quantityDef.type._referencedSection._qualifiedName
  }, [quantityDef])
  const fetchSuggestions = useCallback(input => {
    if (fetchedSuggestionsFor.current === input) {
      return // We've already fetched suggestions for this search string
    }
    // Fetch suggestions
    fetchedSuggestionsFor.current = input
    const query = {}
    if (input !== '') {
      query['entry_name.prefix'] = input
    }
    api.post('entries/query', {
      'owner': 'visible',
      'query': {
        'sections': referencedSectionQualifiedName,
        ...query
      },
      'required': {
        'include': [
          'entry_name',
          'upload_id',
          'entry_id'
        ]
      }
    }, {
      noLoading: true
    }).then(response => {
      setSuggestions(response.data)
    }).catch(raiseError)
  }, [api, raiseError, fetchedSuggestionsFor, setSuggestions, referencedSectionQualifiedName])
  const fetchSuggestionsDebounced = useMemo(() => debounce(fetchSuggestions, 500), [fetchSuggestions])

  useEffect(() => {
    if (!value || value === '') {
      setInputValue('')
      setError(null)
      return
    }
    const resolved = resolveRefAsync(value, archive, context, archive => {
      setEntry(archive.metadata)
      setInputValue(archive.metadata.entry_name)
      setError(null)
    })
    resolved.then(resolved => {
      if (!resolved) {
        setEntry(null)
        setInputValue(value)
        setError('the referenced value does not exist anymore')
      }
    })
  }, [value, api, raiseError, archive, context, setError])

  const getOptionLabel = useCallback(option => option.entry_name, [])
  const getOptionSelected = useCallback((option, value) => option.entry_name === value.entry_name, [])

  const handleValueChange = useCallback((event, value) => {
    if (value?.entry_id && value?.upload_id) {
      value = `../uploads/${value.upload_id}/archive/${value.entry_id}#data`
    } else if (value?.entry_id) {
      value = `../upload/archive/${value.entry_id}#data`
    } else {
      value = undefined
    }
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  const handleInputValueChange = useCallback((event, value) => {
    value = value || event.target.value
    setInputValue(value)
    setError(null)
    fetchSuggestionsDebounced(value)
    if (value === '' && onChange) {
      onChange(undefined)
    }
  }, [setInputValue, fetchSuggestionsDebounced, onChange])

  const handleFocus = useCallback((e) => {
    if (!inputValue) {
      fetchSuggestionsDebounced(inputValue)
    }
  }, [fetchSuggestionsDebounced, inputValue])

  const itemKey = useMemo(() => {
    if (!isNaN(index)) {
      return `${quantityDef.name}:${index}`
    } else {
      return quantityDef.name
    }
  }, [quantityDef, index])
  const {helpDescription, ...otherProps} = getFieldProps(quantityDef)
  return <Autocomplete
    options={suggestions}
    onInputChange={handleInputValueChange}
    onChange={handleValueChange}
    onFocus={handleFocus}
    inputValue={inputValue}
    getOptionLabel={getOptionLabel}
    getOptionSelected={getOptionSelected}
    renderInput={params => {
      return (
        <TextField
          {...params}
          fullWidth variant='filled' size='small'
          {...otherProps}
          {...(value && !entry ? {'data-testid': isWaitingForUpdateTestId} : {})}
          error={!!error}
          helperText={error}
          InputProps={{
            ...params.InputProps,
            endAdornment: inputValue !== '' && (
              <div style={{position: 'absolute', right: 12, top: 'calc(50% - 14px)'}}>
                <ItemButton size="small" itemKey={itemKey} />
              </div>
            )
          }}
        />
      )
    }}
  />
})
ReferenceEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func,
  index: PropTypes.number // additional index used for navigation of higher shaped references
}

export default ReferenceEditQuantity
