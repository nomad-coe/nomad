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
import React, {useCallback, useEffect, useMemo, useState} from 'react'
import PropTypes from 'prop-types'
import {useApi} from '../api'
import {Autocomplete} from '@material-ui/lab'
import {useGlobalMetainfo} from '../archive/metainfo'
import {useDataStore} from '../DataStore'
import {getSchemaInfo, getSectionsInfo} from './SectionSelectDialog'
import {ListItemText} from '@material-ui/core'
import {useSuggestions} from "../../hooks"
import {useErrors} from "../errors"
import { debounce } from 'lodash'

const getOptionSelected = (option, value) => {
  return value ? option.entry_id === value.entry_id && option.value === value.value : false
}

const suggestionQuantities = [{'name': 'entry_name', 'size': 10}, {'name': 'mainfile', 'size': 10}]
const quantitiesAllSet = new Set(['entry_name', 'mainfile'])

const entryFilterOptions = (option) => {
    return `${option.mainfile}|${option.entry_name}|${option.upload_id}|${option.entry_id}`
}

const pathFilterOptions = (option) => {
  return `${option.mainfile}#${option.value}|${option.mainfile}#/${option.value}|${option.mainfile}#${option.label}`
}

const filterOptions = (options, state) => {
  const inputValue = state?.inputValue
  if (inputValue) {
    const parts = inputValue.split("#")
    const entryName = parts?.[0] || ''
    const path = parts?.[1] || ''
    const shownEntries = []
    return options?.filter(option => {
      if (option.menuType === 'entry') {
        const found = entryFilterOptions(option).match(new RegExp(entryName.toLowerCase(), "gi"))
        if (option?.entry_id && found) {
          shownEntries.push(option.entry_id)
        }
        return found
      } else {
        if (option?.entry_id && shownEntries.includes(option.entry_id)) {
          return pathFilterOptions(option).match(new RegExp(path.toLowerCase(), "gi"))
        } else {
          return false
        }
      }
    })
  } else {
    return []
  }
}

function SectionSelectAutocomplete(props) {
  const {renderInput, onValueChanged, value, filtersLocked, onError} = props
  const {api} = useApi()
  const {raiseError} = useErrors()
  const globalMetainfo = useGlobalMetainfo()
  const dataStore = useDataStore()
  const [suggestionInput, setSuggestionInput] = useState('')
  const [entries, setEntries] = useState([])
  const [selectedEntry, setSelectedEntry] = useState(undefined)
  const [suggestions] = useSuggestions(suggestionQuantities, quantitiesAllSet, suggestionInput)
  const [inputValue, setInputValue] = useState('')
  const [options, setOptions] = useState([])
  const [internalError, setInternalError] = useState('')

  useEffect(() => {
    setInternalError('')
    const entryNames = suggestions.filter(suggestion => suggestion.category === 'entry_name').map(suggestion => suggestion.value)
    const mainfiles = suggestions.filter(suggestion => suggestion.category === 'mainfile').map(suggestion => suggestion.value)
    const queries = []
    if (entryNames?.length > 0) {
      queries.push({'entry_name:any': entryNames})
    }
    if (mainfiles?.length > 0) {
      queries.push({'mainfile:any': mainfiles})
    }
    if (value?.entry_id) {
      queries.push({'entry_id': value.entry_id})
    }
    if (queries.length === 0) return
    const requestBody = {
      owner: 'visible',
      exclude: ['quantities', 'sections', 'files'],
      query: {
        and: [
          filtersLocked,
          {or: queries}
        ]
      }
    }
    api.post(`entries/query`, requestBody).then(response => {
      const data = response?.data?.map(entry => {
        return {
          upload_id: entry.upload_id,
          entry_id: entry.entry_id,
          mainfile: entry.mainfile,
          entry_name: entry.entry_name,
          shownValue: '',
          menuType: 'entry'}
      })
      const notFoundEntries = new Set()
      entryNames.forEach(entryName => {
        if (!data.map(entry => entry.entry_name).includes(entryName)) {
          notFoundEntries.add(entryName)
        }
      })
      mainfiles.forEach(mainfile => {
        if (!data.map(entry => entry.mainfile).includes(mainfile)) {
          notFoundEntries.add(mainfile)
        }
      })
      setEntries(data.concat([...notFoundEntries].map(entryName => {
        return {mainfile: entryName, entry_name: entryName, menuType: 'entry'}
      })))
    }).catch(error => {
      raiseError(error)
    })
  }, [api, filtersLocked, raiseError, suggestions, value])

  const handleInputValueChange = useCallback((value) => {
    setSuggestionInput(value.toLowerCase())
  }, [])

  const debouncedHandleInputChange = useMemo(() => {
    return debounce(handleInputValueChange, 700)
  }, [handleInputValueChange])

  const handleInputChange = useCallback((event, value) => {
    if (event?.type === 'change') {
      setInputValue(value)
      debouncedHandleInputChange(value)
    }
  }, [debouncedHandleInputChange])

  useEffect(() => {
    const newOptions = []
    entries.forEach(entry => {
      newOptions.push(entry)
      if (selectedEntry && entry.entry_id === selectedEntry.entryId) {
        selectedEntry.sections.forEach(section => {
          newOptions.push(section)
        })
      }
    })
    if (options.length !== newOptions.length || !options.every((option, index) => option === newOptions[index])) {
      setOptions(newOptions)
    }
  }, [entries, options, selectedEntry])

  const handleChange = useCallback((_, value) => {
    setInputValue(value?.shownValue || '')
    onValueChanged(value)
  }, [onValueChanged])

  const setSection = useCallback((sections, value, updateInputValue) => {
    setSelectedEntry({entryId: value.entry_id, sections: sections})
    if (updateInputValue) {
      setInputValue(`${value.mainfile}#`)
    } else {
      const section = sections.find(section => section.value === value?.value)
      if (!section) {
        onError('The provided path does not exist')
        setInputValue(value?.archive?.mainfile ? (value?.value && value.value !== '/data' && value.value !== 'data' ? `${value.archive.mainfile}#${value?.value}` : value.archive.mainfile) : '')
      } else {
        setInputValue(section.shownValue)
      }
    }
  }, [onError])

  const loadSections = useCallback((value, updateInputValue = false) => {
    if (filtersLocked.entry_type?.includes('Schema')) {
      getSchemaInfo(globalMetainfo, value.entry_id)
        .then(sections => setSection(sections, value, updateInputValue))
    } else {
      getSectionsInfo(api, dataStore, filtersLocked['section_defs.definition_qualified_name'], value?.entry_id)
        .then(sections => setSection(sections, value, updateInputValue))
    }
  }, [api, dataStore, filtersLocked, globalMetainfo, setSection])

  useEffect(() => {
    if (!internalError) {
      if (value?.entry_id && value?.value) {
        loadSections(value)
      }
    }
  }, [internalError, loadSections, onError, value])

  const handleEntryClicked = useCallback((event, option) => {
    event.stopPropagation()
    if (selectedEntry !== undefined && option.entry_id === selectedEntry.entryId) {
      setSelectedEntry(undefined)
    } else {
      loadSections(option, true)
    }
  }, [loadSections, selectedEntry])

  if (options === undefined) {
    return null
  }

  const renderOption = (option, props) => {
    if (option.menuType === 'entry') {
      if (option?.entry_id) {
        return <ListItemText
          key={option.entry_id}
          primary={option?.entry_name || option.mainfile}
          secondary={`upload id: ${option.upload_id}`}
          onClick={(event) => handleEntryClicked(event, option)}
          data-testid={'section-select-entry-activated'}/>
      } else {
        return <ListItemText
          key={option?.entry_name || option.mainfile}
          secondary={option?.entry_name || option.mainfile}
          onClick={(event) => event.stopPropagation()}
          data-testid={'section-select-entry-deactivate'}/>
      }
    } else {
      return <ListItemText
        key={`${option?.entry_id}:${option?.value}`}
        inset={true}
        primary={option?.label}
        data-testid={'section-select-path'}/>
    }
  }

  return <Autocomplete
      value={value || null}
      getOptionSelected={getOptionSelected}
      options={options}
      onInputChange={handleInputChange}
      onChange={handleChange}
      getOptionLabel={(option) => option.shownValue || ''}
      filterOptions={filterOptions}
      renderInput={renderInput}
      renderOption={renderOption}
      inputValue={inputValue}
      freeSolo
  />
}
SectionSelectAutocomplete.propTypes = {
  // The standard autocomplete renderInput
  renderInput: PropTypes.func,
  // The event when the selected section is changed
  onValueChanged: PropTypes.func,
  // The default selected section
  value: PropTypes.object,
  // Determines the fixed search filters to find the corresponding references or the desired sections.
  filtersLocked: PropTypes.object,
  // Fires and passes the error message when an error occurs during fetching the data or resolving the value
  onError: PropTypes.func
}

export default SectionSelectAutocomplete
