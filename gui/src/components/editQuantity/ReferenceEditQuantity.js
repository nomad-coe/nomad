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
import { useApi } from '../api'
import { useErrors } from '../errors'
import {
  Box, Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton, makeStyles,
  TextField,
  Tooltip
} from '@material-ui/core'
import { useEntryStore } from '../entry/EntryContext'
import {ItemButton, useLane} from '../archive/Browser'
import { getFieldProps } from './StringEditQuantity'
import { refType, resolveNomadUrl, getDisplayLabel } from '../../utils'
import AddIcon from '@material-ui/icons/AddCircle'
import { getUrlFromDefinition, QuantityMDef } from '../archive/metainfo'
import EditIcon from '@material-ui/icons/Edit'
import SectionSelectDialog from '../uploads/SectionSelectDialog'
import DialogContentText from '@material-ui/core/DialogContentText'
import AutoComplete from '@material-ui/lab/Autocomplete'
import {useDataStore} from '../DataStore'
import SectionSelectAutocomplete from '../uploads/SectionSelectAutocomplete'
import {Link} from "react-router-dom"
import DetailsIcon from '@material-ui/icons/ArrowForward'
import OverwriteExistingFileDialog from './OverwriteExistingFileDialog'
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

const referenceEditQuantityContext = React.createContext(undefined)

function useReferenceEditQuantityContext() {
  return useContext(referenceEditQuantityContext)
}

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%',
    minWidth: 400
  }
}))

function useReferecedSectionDef(quantityDef) {
  return useMemo(() => {
    const referencedDefinition = quantityDef.type._referencedDefinition
    return referencedDefinition.m_def === QuantityMDef ? referencedDefinition._section : referencedDefinition
  }, [quantityDef])
}

const CreateNewReferenceDialog = React.memo(({allEntryDataSections, open, onSuccess, onFailed, onCanceled}) => {
  const classes = useStyles()
  const {deploymentUrl, uploadId} = useEntryStore('*') || {}
  const {user, api} = useApi()
  const {raiseError} = useErrors()
  const [value, setValue] = useState('')
  const [suggestions, setSuggestions] = useState([])
  const [selectedUpload, setSelectedUpload] = useState(null)
  const [askForOverwrite, setAskForOverwrite] = useState(false)
  const [selectedSection, setSelectedSection] = useState()

  useEffect(() => {
    setSelectedSection(allEntryDataSections?.[0])
    const prepareUploads = async () => {
      const response = await api.get(`/uploads?is_published=false&page_size=10000`)
      const uploads = response.data.map(upload => {
        const uploadName = uploadId === upload.upload_id
          ? (upload.upload_name ? `${upload.upload_name} (This upload)` : 'This upload')
          : (upload.upload_name ? upload.upload_name : `upload-id: ${upload.upload_id}`)
        if (user.sub === upload.main_author) {
          return {label: uploadName, upload_id: upload.upload_id, main_author: upload.main_author, group: 'My uploads'}
        } else {
          return {label: uploadName, upload_id: upload.upload_id, main_author: upload.main_author}
        }
      })
      const myUploads = uploads.filter(upload => upload.group !== undefined)
      const othersUploads = uploads.filter(upload => upload.group === undefined)
      if (othersUploads.length > 0) {
        const authorsID = new Set(othersUploads.map(upload => upload.main_author))
        const response = await api.get(`users?user_id=${Array.from(authorsID).join('&user_id=')}`)
        const users = response['data']
        const usersInfo = {}
        users.forEach(user => (usersInfo[user.user_id] = {name: user.name, affiliation: user.affiliation}))
        othersUploads.forEach(upload => {
          upload['group'] = `Owned by ${(usersInfo[upload.main_author]?.affiliation ? `${usersInfo[upload.main_author]?.name} (${usersInfo[upload.main_author]?.affiliation})` : usersInfo[upload.main_author]?.name)}`
        })
        setSuggestions(myUploads.concat(othersUploads))
      } else {
        setSuggestions(uploads)
      }
    }
    if (user?.sub && open) {
      prepareUploads()
    }
  }, [api, raiseError, uploadId, user?.sub, open, allEntryDataSections])

  const createNewEntry = useCallback((fileName, overwrite = false) => {
    const archive = {
      data: {
        m_def: getUrlFromDefinition(selectedSection, {deploymentUrl, uploadId: selectedUpload.upload_id}, true)
      }
    }
    return new Promise((resolve, reject) => {
      api.put(`uploads/${selectedUpload.upload_id}/raw/?file_name=${fileName}.archive.json&wait_for_processing=true&overwrite_if_exists=${overwrite}`, archive)
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
          reject(error)
        })
    })
  }, [selectedSection, deploymentUrl, selectedUpload?.upload_id, api])

  const handleOverwrite = useCallback((data) => {
    createNewEntry(data, true)
      .then(response => {
        onSuccess(response.processing)
      })
      .catch(error => {
        onFailed(new Error(error))
      })
    setAskForOverwrite(false)
  }, [createNewEntry, onFailed, onSuccess])

  const handleCreateClick = useCallback(() => {
    createNewEntry(value)
      .then(response => {
        onSuccess(response.processing)
      })
      .catch(error => {
        if (error.apiMessage === "The provided path already exists and overwrite_if_exists is set to False.") {
          setAskForOverwrite(true)
        } else {
          onFailed(new Error(error))
        }
      })
  }, [createNewEntry, onFailed, onSuccess, value])

  useEffect(() => {
    if (suggestions.length > 0) {
      const thisUpload = suggestions.find(upload => uploadId === upload.upload_id)
      setSelectedUpload(thisUpload)
    }
  }, [suggestions, uploadId])

  return <React.Fragment>
    <Dialog classes={{paper: classes.dialog}} open={open} disableEscapeKeyDown data-testid='create-reference-dialog'>
      <DialogTitle>{selectedSection?.name ? `Create new reference of type ${selectedSection.name}` : `Create new reference`}</DialogTitle>
      <DialogContent>
        <AutoComplete
          options={suggestions}
          style={{width: '100%', paddingBottom: 10}}
          onChange={(event, value) => setSelectedUpload(value)}
          value={selectedUpload}
          getOptionSelected={(option, value) => option.upload_id === value.upload_id}
          groupBy={(option) => option.group}
          getOptionLabel={(option) => option.label}
          renderInput={(params) => <TextField {...params} label="Target upload" variant='filled' />}
        />
        {allEntryDataSections?.length > 1 ? <AutoComplete
          options={allEntryDataSections}
          style={{width: '100%', paddingBottom: 10}}
          onChange={(event, value) => setSelectedSection(value)}
          value={selectedSection}
          getOptionSelected={(option, value) => option._qualifiedName === value._qualifiedName}
          getOptionLabel={(option) => option.name}
          renderInput={(params) => <TextField {...params} label="Inheriting classes" variant='filled' />}
        /> : ''}
        <TextField
          style={{width: '100%'}}
          label={'Name'}
          value={value}
          variant='filled'
          onChange={event => setValue(event.target.value)}
          data-testid='new-reference-name'
        />
        <DialogContentText>
          {value ? `File name: ${value}.archive.json` : ''}
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={() => onCanceled()} color="secondary">
          Cancel
        </Button>
        <Button onClick={handleCreateClick} disabled={!value || !selectedUpload?.upload_id} color="secondary">
          Create
        </Button>
      </DialogActions>
    </Dialog>
    <OverwriteExistingFileDialog
      open={askForOverwrite}
      data={value}
      onOverwrite={handleOverwrite}
      onCancel={() => setAskForOverwrite(false)}
    />
  </React.Fragment>
})
CreateNewReferenceDialog.propTypes = {
  allEntryDataSections: PropTypes.object.isRequired,
  open: PropTypes.bool,
  onSuccess: PropTypes.func,
  onFailed: PropTypes.func,
  onCanceled: PropTypes.func
}

export const ItemLink = React.forwardRef(function ItemLink({itemKey, ...props}, ref) {
  const reference = useReferenceEditQuantityContext()
  const path = reference?.value
  return <Link
    {...props}
    data-testid={`item:${itemKey}`}
    to={`/user/uploads/upload/id/${reference?.archive?.upload_id}/entry/id/${reference?.archive?.entry_id}/data/${path}`}
  />
})
ItemLink.propTypes = {
  itemKey: PropTypes.string.isRequired
}

const ReferenceEditQuantity = React.memo(function ReferenceEditQuantity(props) {
  const {url} = useEntryStore('*') || {}
  const {quantityDef, value, onChange, index} = props
  const [entry, setEntry] = useState(null)
  const [open, setOpen] = useState(false)
  const [createEntryDialogOpen, setCreateEntryDialogOpen] = useState(false)
  const {user, api} = useApi()
  const {raiseError} = useErrors()
  const [error, setError] = useState({error: undefined, processingError: undefined})
  const lane = useLane()
  const dataStore = useDataStore()
  const referencedSectionDef = useReferecedSectionDef(quantityDef)
  const referencedSectionQualifiedName = useMemo(() => referencedSectionDef?._qualifiedName, [referencedSectionDef])

  useEffect(() => {
    if (!value || value === '') {
      setError({error: null, processingError: null})
      setEntry(null)
      return
    }
    const resolveValue = async () => {
      if (!url) {
        return
      }
      const resolvedUrl = resolveNomadUrl(value, url)
      if (resolvedUrl.type !== refType.archive && resolvedUrl.type !== refType.metainfo) {
        throw new Error(`Expected archive or metainfo reference, got ${resolvedUrl.type} type for ${value}`)
      }
      let query = {
        entry_id: resolvedUrl.entryId
      }
      if (resolvedUrl?.uploadId) {
        query = {upload_id: resolvedUrl.uploadId, ...query}
      }
      const response = await api.post(`entries/archive/query`, {
        owner: 'visible',
        query: query,
        'required': {
          'metadata': {
            'upload_id': '*',
            'entry_id': '*',
            'mainfile': '*',
            'processing_errors': '*'
          }
        }
      }, {
        noLoading: true
      })
      const data = response.data
      if (data.length > 0) {
        const archive = data[0].archive.metadata
        if (archive.processing_errors.length === 0) {
          setEntry({value: value, archive: archive})
          setError({error: null, processingError: null})
        } else {
          setEntry(null)
          setError({error: null, processingError: 'There are some processing errors in the referenced value'})
        }
      } else {
        setEntry(null)
        setError({error: 'The referenced value does not exist anymore', processingError: null})
      }
    }
    resolveValue()
  }, [api, url, value])

  const getReferencePath = useCallback((value) => {
    if (value?.entry_id && value?.upload_id && value?.path) {
      return `../uploads/${value.upload_id}/archive/${value.entry_id}#${value.path}`
    } else if (value?.entry_id && value?.upload_id) {
      return `../uploads/${value.upload_id}/archive/${value.entry_id}#data`
    } else if (value?.entry_id) {
      return `../upload/archive/${value.entry_id}#data`
    } else {
      return undefined
    }
  }, [])

  const changeValue = useCallback((value) => {
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  const handleValueChange = useCallback((value) => {
    const referencePath = getReferencePath(value?.data)
    changeValue(referencePath)
    setOpen(false)
  }, [changeValue, getReferencePath])

  const itemKey = useMemo(() => {
    if (!isNaN(index)) {
      return `${quantityDef.name}:${index}`
    } else {
      return quantityDef.name
    }
  }, [quantityDef, index])

  const {helpDescription, showSectionLabel, ...otherProps} = getFieldProps(quantityDef)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleSuccess = useCallback((value) => {
    setCreateEntryDialogOpen(false)
    const referencePath = getReferencePath({
      entry_name: value.entry.mainfile,
      upload_id: value.upload_id,
      entry_id: value.entry_id
    })
    changeValue(referencePath)
    setEntry({value: referencePath.split('#')[1], archive: value.entry})
  }, [changeValue, getReferencePath])

  const handleFailed = useCallback((error) => {
    setCreateEntryDialogOpen(false)
    raiseError(error)
  }, [raiseError])

  const handleCanceled = useCallback(() => {
    setCreateEntryDialogOpen(false)
  }, [])

  const filtersLocked = useMemo(() => ({'section_defs.definition_qualified_name': [referencedSectionQualifiedName]}), [referencedSectionQualifiedName])

  const addIconButtonToEndAdornment = useCallback((endAdornment, actions) => {
    const children = React.Children.toArray(endAdornment.props.children)
    actions.forEach(iconButton => {
      children.push(iconButton)
    })
    return React.cloneElement(endAdornment, {}, children)
  }, [])

  const allEntryDataSections = useMemo(() => {
    const inheritingSections = dataStore.getAllInheritingSections(referencedSectionDef)
    const allSections = [referencedSectionDef].concat(inheritingSections)
    return allSections.filter(section => !!section._allBaseSections.find(baseSection => baseSection._qualifiedName === 'nomad.datamodel.data.EntryData'))
  }, [dataStore, referencedSectionDef])

  const actions = useMemo(() => {
    const actions = []
    if (!value && allEntryDataSections?.length > 0) {
      actions.push(<IconButton key={'createAction'} size={'small'} disabled={!user?.sub} onClick={() => setCreateEntryDialogOpen(true)}>
        <Tooltip title="Create and assign a new reference">
          <AddIcon/>
        </Tooltip>
      </IconButton>)
    }
    actions.push(<IconButton key={'editAction'} size={'small'} onClick={() => setOpen(true)}>
      <Tooltip title="Search for the references">
        <EditIcon/>
      </Tooltip>
    </IconButton>)
    if (lane && value && !error.error) {
      actions.push(
        // TODO Disabled this button, because the browser does not correctly update after
        // the navigation. This needs to be fixed first.
        <Box display="none">
          <ItemButton
            key={'navigateToReference'}
            size="small" itemKey={itemKey}
            itemLink={ItemLink}
            icon={<DetailsIcon/>}
          />
        </Box>
      )
      actions.push(<ItemButton key={'navigateAction'} size="small" itemKey={itemKey}/>)
    }
    return actions
  }, [value, lane, itemKey, user, error.error, allEntryDataSections])

  const referencedValue = useMemo(() => {
    const value = entry?.value?.split('#')[1] || ''
    return entry ? {entry_id: entry?.archive?.entry_id, value: value.replace(/^(\/*)(.*)/gi, '$2'), archive: entry?.archive} : null
  }, [entry])

  const handleError = useCallback((error) => {
    setError(previousError => ({...previousError, error: error}))
  }, [])

  if (value && referencedValue === undefined) {
    return null
  }

  return <referenceEditQuantityContext.Provider value={referencedValue}>
    <Box display="flex" flexDirection="row" alignItems="center" >
      <Box flexGrow={1}>
        <SectionSelectAutocomplete
          onValueChanged={handleValueChange}
          value={referencedValue}
          filtersLocked={filtersLocked}
          onError={handleError}
          showSectionLabel={showSectionLabel}
          renderInput={(params) => {
            return (
              <TextField
                {...params}
                variant="filled"
                error={!!(error.error || error.processingError)}
                helperText={error.error || error.processingError}
                InputProps={{
                  ...params.InputProps,
                  endAdornment: addIconButtonToEndAdornment(params.InputProps.endAdornment, actions)
                }}
                {...otherProps}
                label={label}
                placeholder={'search by entry name or file name'}
                data-testid='reference-edit-quantity'
              />
            )
          }}
        />
        <CreateNewReferenceDialog open={createEntryDialogOpen} allEntryDataSections={allEntryDataSections} onSuccess={handleSuccess} onFailed={handleFailed} onCanceled={handleCanceled}/>
      </Box>
      <SectionSelectDialog
        open={open}
        onCancel={() => setOpen(false)}
        onSelectedChanged={handleValueChange}
        selected={entry && {entry_id: entry?.archive?.entry_id, value: entry?.value?.split('#')[1]}}
        filtersLocked={filtersLocked}
      />
    </Box>
  </referenceEditQuantityContext.Provider>
})
ReferenceEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func,
  index: PropTypes.number // additional index used for navigation of higher shaped references
}

export default ReferenceEditQuantity
