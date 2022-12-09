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
import {
  Box, Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton, InputAdornment, makeStyles,
  TextField,
  Tooltip
} from '@material-ui/core'
import { useEntryPageContext } from '../entry/EntryPageContext'
import { ItemButton } from '../archive/Browser'
import { getFieldProps } from './StringEditQuantity'
import { isWaitingForUpdateTestId, refType, resolveNomadUrl } from '../../utils'
import AddIcon from '@material-ui/icons/AddCircle'
import { getUrlFromDefinition, QuantityMDef } from '../archive/metainfo'
import EditIcon from '@material-ui/icons/Edit'
import SectionSelectDialog from '../uploads/SectionSelectDialog'
import DialogContentText from '@material-ui/core/DialogContentText'
import ClearIcon from '@material-ui/icons/Clear'
import AutoComplete from '@material-ui/lab/Autocomplete'

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%',
    minWidth: 400
  }
}))

function getReferencedSection(quantityDef) {
  const referencedDefinition = quantityDef.type._referencedDefinition
  const referencedSection = referencedDefinition.m_def === QuantityMDef ? referencedDefinition._section : referencedDefinition
  return referencedSection
}

const CreateNewReference = React.memo(({quantityDef, onSuccess, onFailed}) => {
  const classes = useStyles()
  const {deploymentUrl, uploadId} = useEntryPageContext('*')
  const [open, setOpen] = useState(false)
  const {user, api} = useApi()
  const {raiseError} = useErrors()
  const [value, setValue] = useState('')
  const [suggestions, setSuggestions] = useState([])
  const [selectedUpload, setSelectedUpload] = useState(null)
  const [askForOverwrite, setAskForOverwrite] = useState(false)

  useEffect(() => {
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
  }, [api, raiseError, uploadId, user?.sub, open])

  const createNewEntry = useCallback((fileName, overwrite = false) => {
    const archive = {
      data: {
        m_def: getUrlFromDefinition(getReferencedSection(quantityDef), {deploymentUrl, uploadId: selectedUpload.upload_id}, true)
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
  }, [quantityDef, deploymentUrl, api, selectedUpload])

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
    setOpen(false)
  }, [createNewEntry, onFailed, onSuccess, value])

  useEffect(() => {
    if (suggestions.length > 0) {
      const thisUpload = suggestions.find(upload => uploadId === upload.upload_id)
      setSelectedUpload(thisUpload)
    }
  }, [suggestions, uploadId])

  const handleOverwriteYesClicked = () => {
    createNewEntry(value, true)
        .then(response => {
          onSuccess(response.processing)
        })
        .catch(error => {
          onFailed(new Error(error))
        })
    setAskForOverwrite(false)
    setOpen(false)
  }

  const handleOverwriteNoClicked = () => {
    setAskForOverwrite(false)
    setOpen(false)
  }

  return <Box>
    <IconButton disabled={!user?.sub} onClick={() => setOpen(true)}>
      <Tooltip title="Create and assign a new reference">
        <AddIcon/>
      </Tooltip>
    </IconButton>
    <Dialog classes={{paper: classes.dialog}} open={open} disableEscapeKeyDown data-testid='create-reference-dialog'>
      <DialogTitle>Create new reference</DialogTitle>
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
        <Button onClick={() => setOpen(false)} color="secondary">
          Cancel
        </Button>
        <Button onClick={handleCreateClick} disabled={!value || !selectedUpload?.upload_id} color="secondary">
          Create
        </Button>
      </DialogActions>
    </Dialog>
    <Dialog classes={{paper: classes.dialog}} open={askForOverwrite} disableEscapeKeyDown data-testid='overwrite-reference-dialog'>
      <DialogTitle>Overwrite file</DialogTitle>
      <DialogContent>
        <DialogContentText>
          {`There is already a file with the same name in this path. Do you want to overwrite the file?`}
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={handleOverwriteNoClicked} color="secondary">
          No
        </Button>
        <Button onClick={handleOverwriteYesClicked} color="secondary">
          Yes
        </Button>
      </DialogActions>
    </Dialog>
  </Box>
})
CreateNewReference.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  onSuccess: PropTypes.func,
  onFailed: PropTypes.func
}

const ReferenceEditQuantity = React.memo(function ReferenceEditQuantity(props) {
  const {archive, url} = useEntryPageContext('*')
  const {quantityDef, value, onChange, index} = props
  const [entry, setEntry] = useState(null)
  const [open, setOpen] = useState(false)
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [inputValue, setInputValue] = useState('')
  const [error, setError] = useState()

  const referencedSectionQualifiedName = useMemo(() => {
    const referencedSection = getReferencedSection(quantityDef)
    return referencedSection._qualifiedName
  }, [quantityDef])

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
        const shownValue = resolvedUrl?.path && resolvedUrl.path !== '/data' ? `${entry.mainfile}#${resolvedUrl.path}` : entry.mainfile
        setInputValue(shownValue)
        setError(null)
      } else {
        setEntry(null)
        setInputValue(archive.metadata.mainfile)
      }
    }
    resolveValue()
      .catch(() => {
        setEntry(null)
        setError('the referenced value does not exist anymore')
      })
  }, [value, api, raiseError, archive, url, setError])

  const changeValue = useCallback((value) => {
    if (value?.entry_id && value?.upload_id && value?.path) {
      value = `../uploads/${value.upload_id}/archive/${value.entry_id}#${value.path}`
    } else if (value?.entry_id && value?.upload_id) {
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

  const handleValueChange = useCallback((value) => {
    changeValue(value?.data)
    setOpen(false)
  }, [changeValue])

  const handleClearReference = useCallback(() => {
    setInputValue(undefined)
    changeValue(undefined)
  }, [changeValue])

  const itemKey = useMemo(() => {
    if (!isNaN(index)) {
      return `${quantityDef.name}:${index}`
    } else {
      return quantityDef.name
    }
  }, [quantityDef, index])

  const {helpDescription, ...otherProps} = getFieldProps(quantityDef)

  const handleSuccess = useCallback((value) => {
    setInputValue(value.entry.mainfile)
    changeValue({
      entry_name: value.entry.mainfile,
      upload_id: value.upload_id,
      entry_id: value.entry_id
    })
  }, [changeValue])

  const handleFailed = useCallback((error) => {
    raiseError(error)
  }, [raiseError])

  const filtersLocked = useMemo(() => ({'section_defs.definition_qualified_name': [referencedSectionQualifiedName]}), [referencedSectionQualifiedName])

  return <Box display="flex" flexDirection="row" alignItems="center" >
    <Box flexGrow={1}>
      <TextField
          fullWidth variant='filled' size='small'
          {...otherProps}
          {...(value && !entry ? {'data-testid': isWaitingForUpdateTestId} : {})}
          error={!!error}
          helperText={error}
          InputProps={{
            endAdornment: <InputAdornment position="end">
              {!!inputValue && <IconButton onClick={handleClearReference}>
                <Tooltip title="Delete the reference">
                  <ClearIcon/>
                </Tooltip>
              </IconButton>}
              {!inputValue && <CreateNewReference quantityDef={quantityDef} onSuccess={handleSuccess} onFailed={handleFailed}/>}
              <IconButton onClick={() => setOpen(true)}>
                <Tooltip title="Search for the references">
                  <EditIcon/>
                </Tooltip>
              </IconButton>
              {!!inputValue && <ItemButton size="small" itemKey={itemKey} />}
            </InputAdornment>,
            readOnly: true
          }}
          value={inputValue}
          data-testid='reference-edit-quantity'
      />
    </Box>
    <SectionSelectDialog
        open={open}
        onCancel={() => setOpen(false)}
        onSelectedChanged={handleValueChange}
        selected={value && {entry_id: entry?.entry_id, value: value.split('#')[1]}}
        filtersLocked={filtersLocked}
    />
  </Box>
})
ReferenceEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func,
  index: PropTypes.number // additional index used for navigation of higher shaped references
}

export default ReferenceEditQuantity
