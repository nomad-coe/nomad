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
import { Autocomplete, createFilterOptions } from '@material-ui/lab'
import { makeStyles, TextField } from '@material-ui/core'
import { useEntryPageContext } from '../entry/EntryPageContext'
import { ItemButton } from '../archive/Browser'
import { getFieldProps } from './StringEditQuantity'
import { isWaitingForUpdateTestId, refType, resolveNomadUrl } from '../../utils'
import AddIcon from '@material-ui/icons/AddCircle'

const filter = createFilterOptions()

const useStyles = makeStyles(theme => ({
  icon: {marginRight: theme.spacing(0.5)}
}))

const ReferenceEditQuantity = React.memo(function ReferenceEditQuantity(props) {
  const styles = useStyles()
  const {uploadId, archive, url} = useEntryPageContext('*')
  const {quantityDef, value, onChange, index} = props
  const [entry, setEntry] = useState(null)
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [inputValue, setInputValue] = useState('')
  const [suggestions, setSuggestions] = useState([])
  const [error, setError] = useState()
  const fetchedSuggestionsFor = useRef()

  const referencedSectionQualifiedNames = useMemo(() => {
    return [...quantityDef.type._referencedSection._allInheritingSections.map(section => section._qualifiedName), quantityDef.type._referencedSection._qualifiedName]
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
    const sections = referencedSectionQualifiedNames?.map(qualifiedName => ({'sections': qualifiedName, ...query}))
    api.post('entries/query', {
      'owner': 'visible',
      'query': {
        'or': sections
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
  }, [api, raiseError, fetchedSuggestionsFor, setSuggestions, referencedSectionQualifiedNames])
  const fetchSuggestionsDebounced = useMemo(() => debounce(fetchSuggestions, 500), [fetchSuggestions])

  useEffect(() => {
    if (!value || value === '') {
      setInputValue('')
      setError(null)
      return
    }
    const resolveValue = async () => {
      const resolvedUrl = resolveNomadUrl(value, url)
      if (resolvedUrl.type !== refType.archive) throw new Error(`Archive reference expected, got ${value}`)
      const response = await api.post(`entries/${resolvedUrl.entryId}/archive/query`, {
        'required': {
          'metadata': {
            'entry_name': '*',
            'upload_id': '*',
            'entry_id': '*',
            'mainfile': '*',
            'processing_errors': '*'
          }
        }
      }, {
        noLoading: true
      })
      const entry = response.data.archive.metadata
      if (entry.processing_errors.length === 0) {
        setEntry(entry)
        setSuggestions([{entry_name: entry.entry_name, upload_id: entry.upload_id, entry_id: entry.entry_id}])
        setInputValue(entry.entry_name)
        setError(null)
      } else {
        setEntry(null)
        setSuggestions([{entry_name: archive.metadata.mainfile, upload_id: entry.upload_id, entry_id: entry.entry_id}])
        setInputValue(archive.metadata.mainfile)
      }
    }
    resolveValue()
      .catch(() => {
        setEntry(null)
        setError('the referenced value does not exist anymore')
      })
  }, [value, api, raiseError, archive, url, setError])

  const getOptionLabel = useCallback(option => option.entry_name, [])
  const getOptionSelected = useCallback((option, value) => {
    if (value?.createNewEntry) {
      return true
    }
    return option.entry_name === value.entry_name
  }, [])

  const changeValue = useCallback((value) => {
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

  const createNewEntry = useCallback((uploadId, fileName) => {
    const archive = {
      data: {
        m_def: quantityDef.type._referencedSection._url || quantityDef.type._referencedSection._qualifiedName
      }
    }
    return new Promise((resolve, reject) => {
      api.put(`uploads/${uploadId}/raw/?file_name=${fileName}.archive.json&wait_for_processing=true`, archive)
        .then(response => {
          // TODO handle processing errors
          if (response?.processing?.entry?.process_status !== 'SUCCESS') {
            let error = 'Failed to create the reference.'
            if (response?.processing?.entry?.errors) {
              error = `${error} Details: ${response?.processing?.entry?.errors}`
            }
            reject(new Error(error))
          } else {
            resolve(response)
          }
        })
        .catch(error => {
          reject(new Error(error))
        })
    })
  }, [api, quantityDef.type._referencedSection._qualifiedName, quantityDef.type._referencedSection._url])

  const handleValueChange = useCallback((event, value) => {
    if (value?.createNewEntry) {
      value.entry_name = `${value.createNewEntry}.archive.json`
      createNewEntry(value.upload_id, value.createNewEntry)
        .then(response => {
          setInputValue(response.processing.entry.mainfile)
          changeValue({
            entry_name: response.processing.entry.mainfile,
            upload_id: response.processing.upload_id,
            entry_id: response.processing.entry_id
          })
        })
        .catch(raiseError)
      return
    }
    changeValue(value)
  }, [changeValue, createNewEntry, raiseError])

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

  const renderOption = useCallback(option => {
    return <>
      {option?.createNewEntry && <AddIcon fontSize="small" color="action" className={styles.icon}/>}
      {option.entry_name}
    </>
  }, [styles])

  const filterOptions = useCallback((options, params) => {
    const filtered = filter(suggestions, params)
    let fileName = params.inputValue
    if (params.inputValue !== '') {
      if (params.inputValue.endsWith('.archive.json')) {
        fileName = params.inputValue
      } else if (params.inputValue.endsWith('.archive')) {
        fileName = `${params.inputValue}.json`
      } else {
        fileName = `${params.inputValue}.archive.json`
      }
    }
    if (params.inputValue !== '' && !filtered.map(suggestion => suggestion['entry_name']).includes(fileName)) {
      filtered.push({
        entry_name: `Create "${fileName}" in the current upload`,
        upload_id: uploadId,
        entry_id: '',
        createNewEntry: params.inputValue
      })
    }
    return filtered
  }, [suggestions, uploadId])

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
    filterOptions={filterOptions}
    renderOption={renderOption}
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
