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
import {
  makeStyles, DialogContent, Dialog, Typography, FormControlLabel,
  Box, MenuItem, ListItemIcon, ListItemText, MenuList, Chip, RadioGroup, Radio, FormControl, Checkbox
} from '@material-ui/core'
import Button from '@material-ui/core/Button'
import DialogActions from '@material-ui/core/DialogActions'
import PropTypes from 'prop-types'
import {SearchContext, useSearchContext} from "../search/SearchContext"
import {ui, searchQuantities, apiBase} from "../../config"
import SearchBar from '../search/SearchBar'
import {useApi} from '../api'
import {useUploadPageContext} from './UploadPageContext'
import {useEntryStore} from '../entry/EntryContext'
import {traverse, useGlobalMetainfo} from '../archive/metainfo'
import { defaultFilterGroups, quantityNameSearch } from '../search/FilterRegistry'
import SearchResults from '../search/SearchResults'
import {useDataStore} from '../DataStore'
import {pluralize, resolveNomadUrlNoThrow} from "../../utils"
import {Check} from '@material-ui/icons'
import { cloneDeep } from 'lodash'
import {getItemLabelKey} from '../archive/ArchiveBrowser'

// Use the same context as in the global entries search, but with free-text
// query enabled
const searchDialogContext = React.createContext()
const context = cloneDeep(ui?.apps?.options?.entries)
context.search_syntaxes.exclude = undefined

const allFilters = new Set(defaultFilterGroups && (Object.keys(context?.filter_menus?.options))
  .map(filter => {
    const group = defaultFilterGroups?.[filter]
    return group ? Array.from(group) : []
  }).flat())

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%',
    minWidth: 1100,
    minHeight: 400
  },
  dialogContent: {
  },
  resultsTable: {
    overflowY: 'scroll',
    maxHeight: 700
  },
  searchBar: {
    display: 'flex',
    flexGrow: 0,
    zIndex: 1,
    backgroundColor: 'rgba(0, 0, 0, 0.09)'
  },
  filters: {
    marginTop: 10
  }
}))

const shownColumns = [
  'entry_name',
  'entry_type',
  'authors',
  'upload_name'
]

export async function getSectionsInfo(api, dataStore, reference, entry_id) {
  let response
  try {
    response = await api.post(`entries/archive/query`, {
      owner: 'visible',
      query: {
        entry_id: entry_id,
        processed: true
      }
    })
  } catch (e) {
    return []
  }
  const dataArchive = response?.data?.[0]?.archive
  if (!dataArchive) return
  const sections = Object.keys(dataArchive).filter((key) => 'm_def' in dataArchive[key])
  const referencedSubSections = []
  if (sections.length > 0) {
    for (const section of sections) {
      const m_def = dataArchive[section].m_def
      const url = dataArchive.metadata.upload_id ? `${apiBase}/uploads/${dataArchive.metadata.upload_id}/archive/${entry_id}` : undefined
      const dataMetainfoDefUrl = resolveNomadUrlNoThrow(m_def, url)
      const sectionDef = await dataStore.getMetainfoDefAsync(dataMetainfoDefUrl)
      traverse(dataArchive[section], sectionDef, section, (section, sectionDef, path) => {
        const ref = reference && [...reference][0]
        if (ref &&
          (sectionDef._qualifiedName === ref || sectionDef._allBaseSections?.map(section => section._qualifiedName).includes(ref))) {
          const itemLabelKey = getItemLabelKey(sectionDef)
          const name = itemLabelKey && section[itemLabelKey]
          referencedSubSections.push({name: name, upload_id: response?.data?.[0]?.upload_id, entry_id: response?.data?.[0]?.entry_id, path: path})
        }
      })
    }
  }
  const references = referencedSubSections || []
  return references.map(reference => {
    const fullPath = reference?.path && reference.path !== '/data' && reference.path !== 'data' ? `${dataArchive?.metadata?.mainfile}#${reference.path}` : dataArchive?.metadata?.mainfile
    return {
      label: reference.name ? `${reference.name} (./${reference?.path})` : `./${reference?.path}`,
      fullPath: fullPath,
      mainfile: dataArchive?.metadata?.mainfile,
      entry_name: dataArchive?.metadata?.entry_name,
      entry_id: reference.entry_id,
      upload_id: reference.upload_id,
      shownValue: reference.name || fullPath,
      value: reference.path,
      data: reference
    }
  })
}

export async function getSchemaInfo(globalMetainfo, entry_id) {
  const customMetainfo = await globalMetainfo.fetchAllCustomMetainfos(true, {entry_id: entry_id})
  const schemaArchive = customMetainfo[0]?._data
  const sectionDefs = schemaArchive?.definitions?.section_definitions
  return sectionDefs.map(sectionDef => {
    const fullPath = `${schemaArchive?.metadata?.mainfile}#${sectionDef.name}`
    return {
      label: sectionDef.name,
      fullPath: fullPath,
      mainfile: schemaArchive?.metadata?.mainfile,
      entry_name: schemaArchive?.metadata?.entry_name,
      entry_id: schemaArchive.metadata.entry_id,
      value: sectionDef.name,
      shownValue: sectionDef.name || fullPath,
      data: {sectionDef: sectionDef, archive: schemaArchive}
    }
  })
}

const Details = React.memo(({data}) => {
  const {api} = useApi()
  const entry_id = data?.entry_id
  const globalMetainfo = useGlobalMetainfo()
  const dataStore = useDataStore()
  const {useFiltersLockedState} = useSearchContext()
  const filtersLocked = useFiltersLockedState(['section_defs.definition_qualified_name'])
  const [sections, setSections] = useState()

  const {selected, onSelectedChanged} = useContext(searchDialogContext)

  useEffect(() => {
    if (entry_id && globalMetainfo) {
      fetch()
    }
    async function fetch() {
      if (filtersLocked['section_defs.definition_qualified_name']?.has('nomad.metainfo.metainfo.Definition')) {
        const schemas = await getSchemaInfo(globalMetainfo, entry_id)
        setSections(schemas)
      } else {
        const references = await getSectionsInfo(api, dataStore, filtersLocked['section_defs.definition_qualified_name'], entry_id)
        setSections(references)
      }
    }
  }, [api, dataStore, entry_id, filtersLocked, globalMetainfo])

  if (!sections) {
    return ''
  }

  return (
    <MenuList dense>
      {sections.map((section, index) => {
        const isSelected = selected && section.entry_id === selected.entry_id && section.value === selected.value
        return (
          <MenuItem key={`${entry_id}-${index + 1}`} value={section.label} onClick={() => onSelectedChanged(section)}>
            {isSelected && <ListItemIcon>
              <Check data-testid='check-icon'/>
            </ListItemIcon>}
            <ListItemText inset={!isSelected}>
              {section.label}
            </ListItemText>
          </MenuItem>
        )
      })}
    </MenuList>
  )
})
Details.propTypes = {
  data: PropTypes.object.isRequired
}

function SearchBox({open, onCancel, onSelectedChanged, selected}) {
  const classes = useStyles()
  const {user} = useApi()
  const [filter, setFilter] = useState('visible')
  const [unpublished, setUnpublished] = useState(true)
  const [published, setPublished] = useState(true)
  const {
    useSetFilter,
    useFilters,
    useFiltersLocked,
    useUpdateFilter,
    filters: filterList
  } = useSearchContext()
  const filtersLocked = useFiltersLocked()
  const filters = useFilters([...allFilters]
    .filter(
      filter => filter !== 'visibility' &&
      filter !== 'processed' &&
      filter !== 'upload_id' &&
      filter !== 'published' &&
      filter !== 'main_author.user_id' &&
      !filtersLocked[filter]
    ))
  const updateFilter = useUpdateFilter()
  const uploadContext = useUploadPageContext()
  const entryContext = useEntryStore()
  const {uploadId} = uploadContext || entryContext
  const setVisibilityFilter = useSetFilter('visibility')
  const setUploadIdFilter = useSetFilter('upload_id')
  const setProcessedFilter = useSetFilter('processed')
  const setPublishedFilter = useSetFilter('published')
  const setMainAuthorFilter = useSetFilter('main_author.user_id')

  const handleCancel = useCallback(() => {
    onCancel()
  }, [onCancel])

  useEffect(() => {
    setProcessedFilter(true)
    setPublishedFilter(published === unpublished ? undefined : published)
    if (filter === 'visible') {
      setVisibilityFilter(undefined)
      setMainAuthorFilter(undefined)
      setUploadIdFilter(undefined)
    } else if (filter === 'mine') {
      setVisibilityFilter('user')
      setMainAuthorFilter(undefined)
      setUploadIdFilter(undefined)
    } else if (filter === 'coauthors') {
      setVisibilityFilter(undefined)
      setMainAuthorFilter([user.sub], {queryMode: 'none'})
      setUploadIdFilter(undefined)
    } else if (filter === 'upload') {
      setVisibilityFilter(undefined)
      setMainAuthorFilter(undefined)
      setUploadIdFilter(uploadId)
    }
  }, [filter, published, unpublished, setUploadIdFilter, setMainAuthorFilter, setVisibilityFilter, setPublishedFilter, setProcessedFilter, uploadId, user])

  const contextValue = useMemo(() => ({
    selected: selected,
    onSelectedChanged: onSelectedChanged
  }), [onSelectedChanged, selected])

  const definedFilters = useMemo(() => {
    return Object.values(filters).filter(value => {
      return value !== undefined
    })
  }, [filters])

  const handleResetSearch = () => {
    return Object.keys(filters).forEach(key => {
      updateFilter([key, undefined])
    })
  }

  // Customize the search bar suggestions priority.
  const quantities = useMemo(() => {
    let quantities = new Set([
      'entry_name',
      'mainfile',
      'text_search_contents',
      'results.material.elements',
      'results.material.chemical_formula_hill',
      'results.material.chemical_formula_anonymous',
      'results.material.structural_type',
      'results.material.symmetry.structure_name',
      'results.method.simulation.program_name',
      'authors.name'
    ])
    const allQuantities = [...filterList]
    allQuantities
      .filter(q => !quantities.has(q) && searchQuantities[q]?.suggestion)
      .forEach(q => quantities.add(q))
    quantities = [...quantities].map((name) => ({name, size: 5}))
    quantities.push({name: quantityNameSearch})
    return quantities
  }, [filterList])

  return <searchDialogContext.Provider value={contextValue}>
    <Dialog classes={{paper: classes.dialog}} open={open} data-testid='section-select-dialog'>
      <DialogContent className={classes.dialogContent}>
        <SearchBar quantities={quantities} className={classes.searchBar} />
        <Typography className={classes.filters} variant="body1">
          Filters
        </Typography>
        <Box display='grid'>
          <FormControl>
            <RadioGroup row>
              <FormControlLabel
                value={'visible'}
                control={(
                  <Radio
                    checked={filter === 'visible'}
                    onClick={() => setFilter('visible')}
                  />
                )}
                label='All'
              />
              <FormControlLabel
                value={'mine'}
                control={(
                  <Radio
                    checked={filter === 'mine'}
                    onClick={() => setFilter('mine')}
                  />
                )}
                label='Only mine'
              />
              <FormControlLabel
                value={'coauthors'}
                control={(
                  <Radio
                    checked={filter === 'coauthors'}
                    onClick={() => setFilter('coauthors')}
                  />
                )}
                label='Only from others'
              />
              <FormControlLabel
                value={'upload'}
                control={(
                  <Radio
                    checked={filter === 'upload'}
                    onClick={() => setFilter('upload')}
                  />
                )}
                label='Only this upload'
              />
            </RadioGroup>
          </FormControl>
          <Box display='flex'>
            <FormControlLabel
              control={<Checkbox
                onChange={event => setPublished(event.target.checked)}
                color='primary'
                checked={published}
                name='Published'
              />}
              label='Published'
            />
            <FormControlLabel
              control={<Checkbox
                onChange={event => setUnpublished(event.target.checked)}
                color='primary'
                checked={unpublished}
                name='Unpublished'
              />}
              label='Unpublished'
            />
          </Box>
          {definedFilters.length > 0 && <Chip label={`and ${definedFilters.length} more ${pluralize('filter', definedFilters.length, false)}`} color="primary" onDelete={() => handleResetSearch()}/>}
        </Box>
        <div className={classes.resultsTable}>
          <SearchResults
            defaultUncollapsedEntryID={selected?.entry_id}
            multiSelect={false}
            noAction
          />
        </div>
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={handleCancel} color="secondary">
          Cancel
        </Button>
      </DialogActions>
    </Dialog>
  </searchDialogContext.Provider>
}
SearchBox.propTypes = {
  open: PropTypes.bool,
  onCancel: PropTypes.func,
  onSelectedChanged: PropTypes.func,
  selected: PropTypes.object
}

function SectionSelectDialog(props) {
  const {open, onSelectedChanged, selected, onCancel, filtersLocked} = props
  const columns = context?.columns
  const rows = context?.rows
  columns.selected = shownColumns
  rows.details = {enabled: true, render: Details}
  rows.actions = {enabled: false}

  if (!open) return null

  return <SearchContext
    resource={context?.resource}
    initialPagination={context?.pagination}
    initialColumns={columns}
    initialRows={rows}
    initialFilterMenus={context?.filter_menus}
    initialFiltersLocked={filtersLocked}
    initialSearchSyntaxes={context?.search_syntaxes}
    id='sectionselect'
  >
    <SearchBox open={open} onCancel={onCancel} onSelectedChanged={onSelectedChanged} selected={selected}/>
  </SearchContext>
}
SectionSelectDialog.propTypes = {
  open: PropTypes.bool,
  onSelectedChanged: PropTypes.func,
  selected: PropTypes.object,
  onCancel: PropTypes.func,
  filtersLocked: PropTypes.object
}

export default SectionSelectDialog
