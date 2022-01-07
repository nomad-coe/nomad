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
  makeStyles, DialogTitle, DialogContent, Dialog, IconButton, Tooltip, Divider,
  Typography, TextField, Box
} from '@material-ui/core'
import DialogContentText from '@material-ui/core/DialogContentText'
import EditIcon from '@material-ui/icons/Edit'
import Button from '@material-ui/core/Button'
import DialogActions from '@material-ui/core/DialogActions'
import {uploadPageContext} from '../uploads/UploadPage'
import PropTypes from 'prop-types'
import AutoComplete from '@material-ui/lab/Autocomplete'
import {useApi} from '../api'
import {useErrors} from '../errors'
import {useLocation} from 'react-router-dom'
import {getUrl} from '../nav/Routes'
import AddIcon from '@material-ui/icons/Add'
import {Datatable, DatatableTable} from '../datatable/Datatable'
import DeleteIcon from '@material-ui/icons/Delete'
import OpenInNewIcon from '@material-ui/icons/OpenInNew'
import Quantity from '../Quantity'
import CloseIcon from '@material-ui/icons/Close'
import CheckIcon from '@material-ui/icons/Check'

export const editMetaDataDialogContext = React.createContext()
export const commentContext = React.createContext()
export const referencesContext = React.createContext()
export const datasetsContext = React.createContext()

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%'
  },
  comment: {
    width: '100%'
  },
  datasets: {
    display: 'flex',
    justifyContent: 'left',
    flexWrap: 'wrap',
    listStyle: 'none',
    padding: 'initial'
  }
}))

function EditComments() {
  const classes = useStyles()
  const {data, setComment, setIsCommentChanged} = useContext(editMetaDataDialogContext)
  const [newComment, setNewComment] = useState('')
  const [defaultComment, setDefaultComment] = useState('')

  useEffect(() => {
    if (data?.data?.length > 0) {
      let _comment = data?.data[0]?.entry_metadata?.comment
      setNewComment(_comment)
      setDefaultComment(_comment)
    }
  }, [data])

  const handleTextFieldChange = (event) => {
    let _newComment = event.target.value
    setNewComment(_newComment)
    setComment(_newComment)
    let isSame = _newComment === defaultComment
    setIsCommentChanged(!isSame)
  }

  return <Box display={'block'}>
    <Box marginBottom={1} marginTop={1}>
      <Typography variant='subtitle2'>
        Comments
      </Typography>
    </Box>
    <TextField className={classes.comment} id='outlined-multiline-static' label='Comment'
      multiline rows={6} value={newComment} variant='filled' size='small'
      onChange={handleTextFieldChange}
    />
  </Box>
}

function EditReferences() {
  const {data, api, setIsReferencesChanged, setReferences} = useContext(editMetaDataDialogContext)
  const [newReference, setNewReference] = useState('')
  const [validation, setValidation] = useState('')
  const [newReferences, setNewReferences] = useState([])
  const [defaultReferences, setDefaultReferences] = useState([])
  const [edit, setEdit] = useState({index: -1, value: '', validation: ''})

  const columns = [
    {key: '',
      align: 'left',
      render: reference => {
        let index = newReferences.indexOf(reference)
        if (edit.index === index) {
          return <TextField error={edit.validation !== ''} style={{width: '100%'}} label='' defaultValue={reference}
            helperText={edit.validation !== '' && edit.validation} variant='filled' size='small'
            inputRef={component => component && component.focus()} onChange={handleEditTextFieldChange(index)}/>
        } else {
          return <Box maxWidth='300px' whiteSpace='nowrap' textOverflow='ellipsis' overflow='hidden'>
            <Quantity quantity='data' noLabel noWrap data={{data: reference}} withClipboard />
          </Box>
        }
      }
    }
  ]

  const checkChanges = useCallback((references) => {
    let isSame = defaultReferences.length === references.length &&
        defaultReferences.every(_reference => references.includes(_reference))
    setIsReferencesChanged(!isSame)
  }, [defaultReferences, setIsReferencesChanged])

  const validateAPI = useCallback((value) => {
    let query = {metadata: {references: value}, verify_only: true}
    let error = ''
    api.post(`uploads/${data.upload.upload_id}/edit`, query)
      .catch(err => {
        error = err.apiMessage[0].msg
      })
    return error
  }, [api, data])

  const validate = useCallback((value, index) => {
    if (index !== '' && newReferences.includes(value) && value !== newReferences[index]) return 'Duplicated reference'
    if (index !== '' && value === '') return 'The reference can not be empty'
    if (index === '' && newReferences.includes(value)) return 'The URL is already in the references'
    if ([' ', '\t', '\n'].some(whitespace => value.includes(whitespace))) return 'Whitespace is not allowed'
    return ''
  }, [newReferences])

  const handleEditTextFieldChange = (index) => (event) => {
    let newValue = event.target.value
    setEdit({index: index, value: newValue, validation: validate(newValue, index)})
  }

  useEffect(() => {
    if (data?.data?.length > 0) {
      let _references = data?.data[0]?.entry_metadata?.references
      if (_references) {
        setNewReferences(_references)
        setDefaultReferences(_references)
      }
    }
  }, [data])

  const handleTextFieldChange = (event) => {
    let _newReference = event.target.value
    setNewReference(_newReference)
    setValidation(validate(_newReference, ''))
  }

  const handleAdd = () => {
    if (newReference) {
      let _validation = validate(newReference, '')
      if (_validation === '') _validation = validateAPI(newReference)
      if (_validation === '') {
        let _newReferences = [...newReferences, newReference]
        setNewReferences(_newReferences)
        setReferences(_newReferences)
        checkChanges(_newReferences)
      } else {
        setValidation(_validation)
      }
    }
  }

  const contextValue = {
    newReferences: newReferences,
    setNewReferences: setNewReferences,
    edit: edit,
    setEdit: setEdit,
    checkChanges: checkChanges
  }

  return <referencesContext.Provider value={contextValue}>
    <Box display={'block'}>
      <Box marginBottom={1} marginTop={3}>
        <Typography variant='subtitle2'>
          References
        </Typography>
      </Box>
      <Box display="flex">
        <TextField style={{width: '100%'}} label='Enter the URL' defaultValue='' onChange={handleTextFieldChange}
          error={validation !== ''} helperText={validation !== '' && validation} variant='filled' size='small'/>
        <IconButton size='small' onClick={handleAdd} disabled={validation !== '' || newReference === ''}>
          <AddIcon fontSize='inherit' />
        </IconButton>
      </Box>
      {columns && newReferences.length > 0 && <Datatable columns={columns} data={newReferences} >
        <DatatableTable actions={ReferencesActions}/>
      </Datatable>}
    </Box>
  </referencesContext.Provider>
}

const ReferencesActions = React.memo((props) => {
  const {data} = props
  const {newReferences, setNewReferences, edit, setEdit, checkChanges} = useContext(referencesContext)
  const {setReferences} = useContext(editMetaDataDialogContext)

  let index = newReferences.indexOf(data)

  const handleRemove = () => {
    const filteredReferences = newReferences.filter(reference => !(reference === data))
    setNewReferences(filteredReferences)
    setReferences(filteredReferences)
    checkChanges(filteredReferences)
  }

  const handleEdit = () => {
    setEdit({index: index, value: data, validation: ''})
  }

  const handleOpenLink = () => {
    window.open(data, '_blank')
  }

  const handleApprove = () => {
    let _newReferences = [...newReferences]
    _newReferences[edit.index] = edit.value
    setNewReferences(_newReferences)
    setReferences(_newReferences)
    checkChanges(_newReferences)
    setEdit({index: -1, value: '', validation: ''})
  }

  const handleCancel = () => {
    setEdit({index: -1, value: '', validation: ''})
  }

  if (edit.index === index) {
    return <Box display={'inline-block'}>
      <IconButton size='small' onClick={handleCancel}>
        <Tooltip title='Cancel editing'>
          <CloseIcon />
        </Tooltip>
      </IconButton>
      <IconButton size='small' onClick={handleApprove} disabled={edit.validation !== ''}>
        <Tooltip title='Approve editing'>
          <CheckIcon />
        </Tooltip>
      </IconButton>
    </Box>
  }

  return <Box display={'inline-block'}>
    <IconButton size='small' onClick={handleEdit}>
      <Tooltip title='Edit the reference'>
        <EditIcon />
      </Tooltip>
    </IconButton>
    <IconButton size='small' onClick={handleOpenLink}>
      <Tooltip title='Open in new tab'>
        <OpenInNewIcon />
      </Tooltip>
    </IconButton>
    <IconButton size='small' onClick={handleRemove}>
      <Tooltip title='Remove the reference'>
        <DeleteIcon />
      </Tooltip>
    </IconButton>
  </Box>
})
ReferencesActions.propTypes = {
  data: PropTypes.object.isRequired
}

function EditDatasets() {
  const {data, api, raiseError, setIsDatasetChanged, setDatasets} = useContext(editMetaDataDialogContext)
  const [suggestions, setSuggestions] = useState([])
  const [validation, setValidation] = useState('')
  const [allDatasets, setAllDatasets] = useState([])
  const [newDataset, setNewDataset] = useState({dataset_id: '', dataset_name: ''})
  const [addDataset, setAddDataset] = useState('')
  const [newDatasets, setNewDatasets] = useState([])
  const [isDuplicated, setIsDuplicated] = useState(false)
  const [defaultDatasets, setDefaultDatasets] = useState([])

  const columns = useMemo(() => ([
    {key: '', align: 'left', render: dataset => dataset.dataset_name}
  ]), [])

  useEffect(() => {
    api.get(`/datasets?page_size=${1000}&page=${1}`)
      .then(datasets => {
        setAllDatasets(datasets?.data)
      })
      .catch(err => {
        setAllDatasets([])
        raiseError(err)
      })
  }, [api, raiseError])

  const checkChanges = useCallback((_datasets) => {
    let isSame = defaultDatasets.length === _datasets.length &&
        defaultDatasets.map(dataset => dataset.dataset_id).every(id => _datasets.map(dataset => dataset.dataset_id).includes(id))
    setIsDatasetChanged(!isSame)
  }, [defaultDatasets, setIsDatasetChanged])

  const handleInputChange = useCallback((event, value) => {
    const query = value.toLowerCase()
    const withQueryInName = (allDatasets ? allDatasets.filter(dataset => dataset.dataset_name.toLowerCase().indexOf(query) !== -1) : [])
    setSuggestions(withQueryInName.slice(0, 5))
  }, [allDatasets])

  useEffect(() => {
    if (data?.data?.length > 0) {
      let _datasets = data?.data[0]?.entry_metadata?.datasets
      if (_datasets) {
        setNewDatasets(_datasets)
        setDefaultDatasets(_datasets)
      }
    }
  }, [data])

  const validate = useCallback((value) => {
    if (allDatasets.map(dataset => dataset.dataset_name).includes(value.dataset_name)) return `There is already a dataset with name ${value.dataset_name}`
    if (newDatasets.map(dataset => dataset.dataset_name).includes(value.dataset_name)) return `There is already a dataset with name ${value.dataset_name}`
    if (value.dataset_name[0] === ' ') return `Invalid name for dataset`
    return ''
  }, [allDatasets, newDatasets])

  const handleAutoCompleteChange = (event, value) => {
    if (value && value?.dataset_id) {
      setAddDataset(value)
      setIsDuplicated(newDatasets.map(dataset => dataset.dataset_id).includes(value.dataset_id))
    } else {
      setAddDataset('')
    }
  }

  const handleCreate = () => {
    if (newDataset) {
      if (!newDatasets.map(dataset => dataset.dataset_id).includes(newDataset.dataset_id)) {
        let _newDatasets = [...newDatasets, newDataset]
        setNewDatasets(_newDatasets)
        setDatasets(_newDatasets)
        checkChanges(_newDatasets)
      } else {
        setValidation(validate(newDataset))
      }
    }
  }

  const handleAdd = () => {
    if (addDataset) {
      if (!newDatasets.map(dataset => dataset.dataset_id).includes(addDataset.dataset_id)) {
        let _newDatasets = [...newDatasets, addDataset]
        setNewDatasets(_newDatasets)
        setDatasets(_newDatasets)
        let isSame = defaultDatasets.length === _newDatasets.length &&
            defaultDatasets.map(dataset => dataset.dataset_id).every(id => _newDatasets.map(dataset => dataset.dataset_id).includes(id))
        setIsDatasetChanged(!isSame)
      } else {
        setIsDuplicated(true)
      }
    }
  }

  const handleTextFieldChange = (event) => {
    let _newDataset = {dataset_id: event.target.value, dataset_name: event.target.value}
    setNewDataset(_newDataset)
    setValidation(validate(_newDataset))
  }

  const contextValue = useMemo(() => ({
    newDatasets: newDatasets,
    setNewDatasets: setNewDatasets,
    setDatasets: setDatasets,
    checkChanges: checkChanges
  }), [newDatasets, setNewDatasets, setDatasets, checkChanges])

  return <datasetsContext.Provider value={contextValue}>
    <Box display={'block'}>
      <Box marginBottom={1} marginTop={3}>
        <Typography variant='subtitle2'>
          Datasets
        </Typography>
      </Box>
      <Box display="flex" marginBottom={1}>
        <TextField style={{width: '100%'}} label='Create a new dataset' placeholder='Dataset name' defaultValue='' onChange={handleTextFieldChange}
          error={validation !== ''} helperText={validation !== '' && validation} variant='filled' size='small'/>
        <IconButton size="small" onClick={handleCreate} disabled={validation !== '' || newDataset.dataset_name === ''}>
          <AddIcon fontSize="inherit" />
        </IconButton>
      </Box>
      <Box display="flex">
        <AutoComplete
          style={{width: '100%'}}
          options={suggestions}
          getOptionLabel={option => option.dataset_name}
          onInputChange={handleInputChange}
          onChange={handleAutoCompleteChange}
          renderInput={params => (
            <TextField
              {...params}
              variant='filled' label='Search the datasets' placeholder='Dataset name' argin='normal' fullWidth size='small'
              error={isDuplicated} helperText={isDuplicated && 'The data is already in the selected dataset'}
            />
          )}
        />
        <IconButton size="small" onClick={handleAdd} disabled={isDuplicated || addDataset === ''}>
          <AddIcon fontSize="inherit" />
        </IconButton>
      </Box>
      {columns && newDatasets.length > 0 && <Datatable columns={columns} data={newDatasets} >
        <DatatableTable actions={DatasetsActions}/>
      </Datatable>}
    </Box>
  </datasetsContext.Provider>
}

const DatasetsActions = React.memo((props) => {
  const {data} = props
  const location = useLocation()
  const {newDatasets, setNewDatasets, setDatasets, checkChanges} = useContext(datasetsContext)

  const handleRemove = () => {
    const filteredDatasets = newDatasets.filter(dataset => !(dataset.dataset_id === data.dataset_id))
    setNewDatasets(filteredDatasets)
    setDatasets(filteredDatasets)
    checkChanges(filteredDatasets)
  }

  const handleOpenLink = () => {
    const path = `dataset/id/${data.dataset_id}`
    window.open(getUrl(path, location), '_blank')
  }

  return <Box display={'inline-block'}>
    <IconButton size='small' onClick={handleOpenLink}>
      <Tooltip title="Open in new tab">
        <OpenInNewIcon />
      </Tooltip>
    </IconButton>
    <IconButton size='small' onClick={handleRemove}>
      <Tooltip title="Remove the dataset">
        <DeleteIcon />
      </Tooltip>
    </IconButton>
  </Box>
})
DatasetsActions.propTypes = {
  data: PropTypes.object.isRequired
}

function EditMetaDataDialog({...props}) {
  const {isIcon, selectedEntries} = props
  const classes = useStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const {upload, setUpload, data} = useContext(uploadPageContext)
  const [open, setOpen] = useState(false)
  const [openConfirmDialog, setOpenConfirmDialog] = useState(false)
  const isProcessing = upload?.process_running
  const nSelected = (selectedEntries.upload_id ? upload?.entries : selectedEntries.entry_id.length)

  const [comment, setComment] = useState('')
  const [references, setReferences] = useState([])
  const [datasets, setDatasets] = useState([])
  const [isCommentChanged, setIsCommentChanged] = useState(false)
  const [isReferencesChanged, setIsReferencesChanged] = useState(false)
  const [isDatasetChanged, setIsDatasetChanged] = useState(false)

  const handleOpenDialog = () => {
    setIsCommentChanged(false)
    setIsReferencesChanged(false)
    setIsDatasetChanged(false)
    setOpen(true)
  }

  const handleDiscardChanges = () => {
    setOpenConfirmDialog(false)
    setOpen(false)
  }

  const createNewDatasets = useCallback(() => {
    let promises = []
    let newDatasets = datasets.filter(_dataset => _dataset.dataset_id === _dataset.dataset_name).map(_dataset => _dataset.dataset_name)
    newDatasets.forEach(dataset_name => {
      promises.push(api.post(`/datasets`, {dataset_name: dataset_name}))
    })
    return Promise.all(promises)
  }, [api, datasets])

  const submitChanges = useCallback((metadata) => {
    let requestBody = {metadata: metadata, verify_only: false, owner: 'user'}
    if (selectedEntries.entry_id) requestBody.query = {entry_id: {any: selectedEntries.entry_id}}
    api.post(`uploads/${data.upload.upload_id}/edit`, requestBody)
      .then(results => {
        setUpload(results.data)
        setOpen(false)
      }).catch(err => {
        raiseError(err)
      })
  }, [api, raiseError, data, setUpload, selectedEntries])

  const handleSubmitChanges = () => {
    if (isCommentChanged || isReferencesChanged || isDatasetChanged) {
      let metadata = {}
      if (isCommentChanged) metadata.comment = comment
      if (isReferencesChanged) metadata.references = references
      if (isDatasetChanged) {
        createNewDatasets().then(newDatasets => {
          metadata.datasets = datasets.filter(_dataset => _dataset.dataset_id !== _dataset.dataset_name)
            .map(_dataset => _dataset.dataset_id)
            .concat(newDatasets.map(_dataset => _dataset.dataset_id))
          submitChanges(metadata)
        })
      } else {
        submitChanges(metadata)
      }
    } else {
      setOpen(false)
    }
  }

  const handleConfirm = () => {
    if (isCommentChanged || isReferencesChanged || isDatasetChanged) {
      setOpenConfirmDialog(true)
    } else {
      setOpen(false)
    }
  }

  const contextValue = useMemo(() => ({
    api: api,
    raiseError: raiseError,
    setIsCommentChanged: setIsCommentChanged,
    setIsReferencesChanged: setIsReferencesChanged,
    setIsDatasetChanged: setIsDatasetChanged,
    upload: upload,
    data: data,
    setComment: setComment,
    setReferences: setReferences,
    setDatasets: setDatasets
  }), [api, raiseError, setIsCommentChanged, setIsReferencesChanged, setIsDatasetChanged, upload, data, setComment, setReferences, setDatasets])

  return <editMetaDataDialogContext.Provider value={contextValue}>
    <React.Fragment>
      {isIcon && <IconButton onClick={handleOpenDialog}>
        <Tooltip title="Manage upload members">
          <EditIcon />
        </Tooltip>
      </IconButton>}
      {!isIcon && <Button onClick={handleOpenDialog} variant='contained' color='primary' disabled={isProcessing}>
        {`Edit metadata of all ${upload?.entries} entries`}
      </Button>}
      {open && <Dialog classes={{paper: classes.dialog}} open={open} on disableBackdropClick disableEscapeKeyDown>
        <DialogTitle>Manage upload meta data</DialogTitle>
        <DialogContent>
          <DialogContentText>
            Upload ID: {upload?.upload_id}
            <br/>
            You can add, remove or edit the meta data for the selected entries.
            <br/>
            {nSelected} of {upload?.entries} {upload?.entries === 1 ? 'entry' : 'entries'} is selected.
          </DialogContentText>
          <Divider/>
          <div>
            <EditComments/>
            <EditReferences/>
            <EditDatasets/>
          </div>
        </DialogContent>
        <DialogActions>
          <span style={{flexGrow: 1}} />
          <Button onClick={handleConfirm} color="secondary">
            Cancel
          </Button>
          <Button onClick={handleSubmitChanges} disabled={!(isCommentChanged || isReferencesChanged || isDatasetChanged)} color="secondary">
            Submit
          </Button>
        </DialogActions>
        <Dialog
          open={openConfirmDialog}
          aria-describedby="alert-dialog-description"
        >
          <DialogContent>
            <DialogContentText id="alert-dialog-description">
              Your changes are not submitted yet. Discard changes?
            </DialogContentText>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setOpenConfirmDialog(false)} autoFocus>Cancel</Button>
            <Button onClick={handleDiscardChanges}>Discard</Button>
          </DialogActions>
        </Dialog>
      </Dialog>}
    </React.Fragment>
  </editMetaDataDialogContext.Provider>
}
EditMetaDataDialog.propTypes = {
  isIcon: PropTypes.bool,
  selectedEntries: PropTypes.arrayOf(PropTypes.object)
}

export default EditMetaDataDialog
