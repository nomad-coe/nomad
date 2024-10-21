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
import {cloneDeep} from 'lodash'
import {
  Chip,
  Dialog, DialogContent, IconButton, TextField, Tooltip
} from "@material-ui/core"
import {FreeformSearchContext, useSearchContext} from "../search/SearchContext"
import {ui} from "../../config"
import DialogActions from "@material-ui/core/DialogActions"
import Button from "@material-ui/core/Button"
import SearchPage from "../search/SearchPage"
import SearchIcon from "@material-ui/icons/Search"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"
import {getDisplayLabel, pluralize} from "../../utils"
import Autocomplete from "@material-ui/lab/Autocomplete"
import {ItemButton, useLane} from "../archive/Browser"
import ClearIcon from "@material-ui/icons/Clear"

const context = cloneDeep(ui?.apps?.options?.eln)
const shownColumns = [
  'entry_name',
  'entry_type',
  'authors',
  'upload_name'
]
const rows = context?.rows
const columns = (context?.columns || [])
  .map((column) => ({...column, selected: shownColumns.includes(column.quantity)}))
rows.details = {enabled: false}
rows.actions = {enabled: false}

function SearchDialog({open, filters, pageSize, onCancel, onQueryChanged}) {
  const {filters: queryFilters, useFilters, useSetFilters, useResults, useApiData} = useSearchContext()
  const setFilters = useSetFilters()
  const filterValues = useFilters(queryFilters)
  const {data, setPagination} = useResults()
  const apiData = useApiData()

  const updateFilters = useCallback(() => {
    const newValue = {}
    for (const key in filters) {
      if (Object.hasOwnProperty.call(filters, key)) {
        newValue[key] = Array.isArray(filters[key]) || filters[key] instanceof Set ? new Set(filters[key]) : filters[key]
      }
    }
    if (setPagination && pageSize) {
      setPagination(old => {
        return {...old, page_size: pageSize}
      })
    }
    setFilters(newValue)
  }, [pageSize, setFilters, filters, setPagination])

  useEffect(() => {
    if (open) {
      updateFilters()
    }
  }, [updateFilters, open])

  const newFilters = useMemo(() => {
    const filters = {...filterValues}
    for (const key in filters) {
      if (key in filters && filters[key] === undefined) {
        delete filters[key]
      }
    }
    return filters
  }, [filterValues])

  const results = useMemo(() => {
    if (!data) {
      return undefined
    }
    return data.map(entry => ({entry_id: entry.entry_id, upload_id: entry.upload_id, mainfile: entry.mainfile}))
  }, [data])

  const handleQueryChanged = useCallback(() => {
    if (onQueryChanged) {
      onQueryChanged(newFilters, apiData, results)
    }
  }, [onQueryChanged, newFilters, apiData, results])

  return <Dialog
    open={open}
    PaperProps={{
      style: {
        maxWidth: '1800px',
        maxHeight: '1200px',
        width: '1800px',
        height: '1200px'
      }
    }}
    data-testid='search-dialog'
  >
    <DialogContent>
      <SearchPage/>
    </DialogContent>
    <DialogActions>
      <span style={{flexGrow: 1}} />
      <Button onClick={() => onCancel()} color="secondary">
        Cancel
      </Button>
      <Button onClick={() => handleQueryChanged()} color="secondary" data-testid='search-dialog-ok'>
        OK
      </Button>
    </DialogActions>
  </Dialog>
}
SearchDialog.propTypes = {
  open: PropTypes.bool,
  filters: PropTypes.object,
  pageSize: PropTypes.number,
  onCancel: PropTypes.func,
  onQueryChanged: PropTypes.func
}

function QueryEditQuantity({quantityDef, onChange, value, storeInArchive, index, maxData}) {
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)
  const [open, setOpen] = useState(false)
  const lane = useLane()

  const filters = useMemo(() => value?.filters || {}, [value])

  const handleCancel = useCallback(() => {
    setOpen(false)
  }, [])

  const handleQueryChanged = useCallback((filters, query, results) => {
    if (onChange) {
      const newFilters = {}
      for (const key in filters) {
        if (Object.hasOwnProperty.call(filters, key)) {
          newFilters[key] = Array.isArray(filters[key]) || filters[key] instanceof Set ? [...filters[key]] : filters[key]
        }
      }
      const newValue = {
        owner: query.response.owner,
        query: query.response.query,
        pagination: query.response.pagination,
        filters: newFilters
      }
      if (storeInArchive) {
        newValue.data = results
      }
      onChange(newValue)
    }
    setOpen(false)
  }, [onChange, storeInArchive])

  const tags = useMemo(() => {
    let tags = []
    for (const key in filters) {
      const filterValue = filters[key]
      if (filterValue) {
        if (Array.isArray(filterValue[key])) {
          tags = tags.concat([...filterValue].map(value => ({key: key, value: value, tag: `${key}:${value}`})))
        } else {
          tags = tags.concat({key: key, value: filterValue, tag: `${key}:${filterValue}`})
        }
      }
    }
    return tags
  }, [filters])

  const itemKey = useMemo(() => {
    if (!isNaN(index)) {
      return `${quantityDef.name}:${index}`
    } else {
      return quantityDef.name
    }
  }, [quantityDef, index])

  const handleClearResults = useCallback(() => {
    if (onChange) {
      onChange(undefined)
    }
  }, [onChange])

  const actions = useMemo(() => {
    const actions = []
    if (value) {
      actions.push(
        <IconButton
          key={'clear'}
          color="seconadry"
          size="small"
          onClick={handleClearResults}
        >
          <Tooltip title="Clear results">
            <ClearIcon/>
          </Tooltip>
        </IconButton>
      )
    }
    actions.push(<IconButton
      key={'search'}
      color="seconadry"
      size="small"
      onClick={() => setOpen(true)}
    >
      <Tooltip title="Search dialog">
        <SearchIcon/>
      </Tooltip>
    </IconButton>)
    if (lane && value) {
      actions.push(<ItemButton key={'navigate'} size="small" itemKey={itemKey}/>)
    }
    return actions
  }, [value, lane, handleClearResults, itemKey])

  return <React.Fragment>
    <Autocomplete
      multiple
      open={false}
      options={[]}
      getOptionLabel={(option) => option.tag}
      limitTags={1}
      value={tags}
      inputValue={''}
      renderTags={(value, getTagProps) => {
        return value.map((option, index) => {
          return <Chip
            key={index}
            {...getTagProps({ index })}
            label={option.tag}
            size="small"
            color="primary"
            onDelete={undefined}
          />
        })
      }}
      renderInput={(params) => (
        <TextField
          {...params}
          label={label}
          variant="filled"
          placeholder={value?.results && Array.isArray(value.results) && value.results.length > 0
            ? pluralize('result', value.results.length, true)
            : "no results"}
          InputProps={{
            ...params.InputProps,
            endAdornment: React.cloneElement(params.InputProps.endAdornment, {}, actions)
          }}
        />
      )}
    />
    <FreeformSearchContext
      resource={context?.resource}
      initialPagination={context?.pagination}
      initialColumns={columns}
      initialRows={rows}
      initialMenu={context?.menu}
      initialFiltersLocked={undefined}
      initialFilterValues={filters}
      initialSearchSyntaxes={context?.search_syntaxes}
      id={`queryeditquantity-${quantityDef._qualifiedName}`}
    >
      <SearchDialog
        open={open}
        filters={filters}
        onCancel={handleCancel}
        onQueryChanged={handleQueryChanged}
        pageSize={maxData || 100}
      />
    </FreeformSearchContext>
  </React.Fragment>
}
QueryEditQuantity.propTypes = {
  // The quantity definition
  quantityDef: PropTypes.object,
  // The event when the searched results have been changed
  onChange: PropTypes.func,
  // The searched value
  value: PropTypes.string,
  // To store the search results in the value
  storeInArchive: PropTypes.bool,
  // The index of the quantity which repeats
  index: PropTypes.number,
  // The maximum number of searched data
  maxData: PropTypes.number
}

export default QueryEditQuantity
