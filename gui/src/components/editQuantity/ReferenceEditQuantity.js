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
import React, { useCallback, useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { debounce } from 'lodash'
import { Autocomplete } from '@material-ui/lab'
import { TextField } from '@material-ui/core'
import { useEntryContext } from '../entry/EntryContext'
import { resolveRefAsync } from '../archive/metainfo'
import { ItemButton, useLane } from '../archive/Browser'
import { getFieldProps } from './StringEditQuantity'

const ReferenceEditQuantity = React.memo(function ReferenceEditQuantity(props) {
  const {uploadId} = useEntryContext()
  const {archive} = useEntryContext()
  const {quantityDef, value, onChange, index} = props
  const [entry, setEntry] = useState(null)
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [inputValue, setInputValue] = useState('')
  const [suggestions, setSuggestions] = useState([])
  const {adaptor: {context}} = useLane()
  const referencedSectionQualifiedName = useMemo(() => {
    return quantityDef.type._referencedSection._qualifiedName
  }, [quantityDef])
  const fetchSuggestions = useCallback(input => {
    const query = {
      'upload_id': uploadId
    }
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
      const suggestions = response.data
      setSuggestions(suggestions)
    }).catch(raiseError)
  }, [api, raiseError, setSuggestions, referencedSectionQualifiedName, uploadId])

  useEffect(() => {
    if (!value || value === '') {
      return
    }
    resolveRefAsync(value, archive, context, archive => {
      setEntry(archive.metadata)
      setInputValue(archive.metadata.entry_name)
    })
  }, [value, api, raiseError, archive, context])

  const fetchSuggestionsDebounced = useCallback(
    debounce(fetchSuggestions, 150),
    [fetchSuggestions]
  )

  useEffect(() => {
    fetchSuggestionsDebounced(inputValue)
  }, [fetchSuggestionsDebounced, inputValue])

  const handleValueChange = useCallback((event, value) => {
    if (value?.entry_id) {
      value = `../upload/archive/${value.entry_id}#data`
    } else {
      value = undefined
    }
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  const getOptionLabel = useCallback(option => option.entry_name, [])

  const handleInputValueChange = useCallback((event, value) => {
    value = value || event.target.value
    setInputValue(value)
    if (value === '' && onChange) {
      onChange(undefined)
    }
  }, [setInputValue, onChange])

  const handleFocus = useCallback(() => {
    if (inputValue === '') {
      fetchSuggestionsDebounced('')
    }
  }, [fetchSuggestionsDebounced, inputValue])

  const itemKey = useMemo(() => {
    if (!isNaN(index)) {
      return `${quantityDef.name}:${index}`
    } else {
      return quantityDef.name
    }
  }, [quantityDef, index])

  return <Autocomplete
    disabled={value && !entry}
    options={suggestions}
    onInputChange={handleInputValueChange}
    onChange={handleValueChange}
    onFocus={handleFocus}
    inputValue={inputValue}
    getOptionLabel={getOptionLabel}
    renderInput={params => {
      return (
        <TextField
          {...params}
          fullWidth variant='filled' size='small'
          {...getFieldProps(quantityDef)}
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
