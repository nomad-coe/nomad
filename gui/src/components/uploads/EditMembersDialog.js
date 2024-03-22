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
import {
  Box,
  Checkbox,
  CircularProgress,
  Dialog,
  DialogContent,
  DialogTitle,
  Divider,
  FormControl,
  FormControlLabel,
  IconButton,
  InputLabel,
  MenuItem, Select,
  TextField,
  Tooltip,
  Typography,
  makeStyles
} from '@material-ui/core'
import Button from '@material-ui/core/Button'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContentText from '@material-ui/core/DialogContentText'
import DeleteIcon from '@material-ui/icons/Delete'
import MembersIcon from '@material-ui/icons/People'
import AutoComplete from '@material-ui/lab/Autocomplete'
import { debounce } from 'lodash'
import PropTypes from 'prop-types'
import React, { useCallback, useContext, useMemo, useReducer, useState } from 'react'
import { useApi } from '../api'
import { Datatable, DatatableTable } from '../datatable/Datatable'
import { useErrors } from '../errors'
import { useUploadPageContext } from './UploadPageContext'

export const editMembersDialogContext = React.createContext()

const useInviteUserDialogStyles = makeStyles(theme => ({
  button: {
    marginLeft: theme.spacing(1)
  },
  dialog: {
    width: '100%'
  },
  submitWrapper: {
    margin: theme.spacing(1),
    position: 'relative'
  },
  submitProgress: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12
  }
}))

const InviteUserDialog = React.memo(function InviteUserDialog(props) {
  const classes = useInviteUserDialogStyles()
  const [open, setOpen] = useState(false)
  const [submitting, setSubmitting] = useState(false)
  const [canSubmit, setCanSubmit] = useState(false)
  const [error, setError] = useState(null)
  const [data, setData] = useState({
    first_name: '',
    last_name: '',
    email: '',
    affiliation: ''
  })
  const {api} = useApi()

  const handleClose = useCallback((event, reason) => {
    if (reason !== 'backdropClick') {
      setOpen(false)
    }
  }, [setOpen])

  const handleSubmit = useCallback(() => {
    setSubmitting(true)

    api.inviteUser(data).then(() => {
      setSubmitting(false)
      setOpen(false)
    }).catch(error => {
      let message = '' + error
      try {
        message = JSON.parse(error.request.responseText).detail || '' + error
      } catch (e) {}
      setError(message)
      setSubmitting(false)
      setCanSubmit(false)
    })
  }, [data, setSubmitting, setOpen, setError, api])

  const handleChange = useCallback((key, value) => {
    const valid = value && !Object.keys(data).find(dataKey => !(key === dataKey || data[dataKey]))
    setData({...data, [key]: value})
    setCanSubmit(valid)
  }, [setData, data, setCanSubmit])

  const handleOpen = useCallback(() => {
    setOpen(true)
  }, [setOpen])

  const input = (key, label) => <TextField
    variant="filled"
    label={label}
    value={data[key]}
    onChange={event => handleChange(key, event.target.value)}
    margin="normal"
    fullWidth
  />
  return <React.Fragment>
    <Button className={classes.button}
      onClick={handleOpen}
      color="secondary" disabled={submitting}
    >
      Invite new user
    </Button>
    <Dialog
      classes={{paper: classes.dialog}}
      open={open}
      onClose={handleClose} disableEscapeKeyDown>
      <DialogTitle>Invite a new user to NOMAD</DialogTitle>
      <DialogContent>
        <DialogContentText>
          If you want to add a user as co-author or share your data with someone that
          is not already a NOMAD user, you can invite this person here. We need just a few
          details about this person. After your invite, the new user will receive an
          Email that allows her to set a password and further details. Anyhow, you will
          be able to add the user as co-author or someone to share with immediately after the
          invite.
        </DialogContentText>
        {error && <DialogContentText color="error">
          {error}
        </DialogContentText>}
        {input('email', 'Email')}
        {input('first_name', 'First name')}
        {input('last_name', 'Last name')}
        {input('affiliation', 'Affiliation')}
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} disabled={submitting}>
          Cancel
        </Button>
        <div className={classes.submitWrapper}>
          <Button onClick={handleSubmit} color="primary" disabled={!canSubmit}>
            Submit
          </Button>
          {submitting && <CircularProgress size={24} className={classes.submitProgress} />}
        </div>
      </DialogActions>
    </Dialog>
  </React.Fragment>
})
InviteUserDialog.propTypes = {}

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%'
  }
}))

function MembersTable() {
  const {members, setIsChanged} = useContext(editMembersDialogContext)
  const forceUpdate = useReducer(bool => !bool)[1]

  const columns = [
    {key: 'Name', align: 'left', render: member => member.name},
    {key: 'Affiliation', align: 'left', render: member => member.affiliation},
    {
      key: 'Role',
      align: 'left',
      render: member => (member.role === 'Main author' ? member.role
        : <Select value={member.role}
          onChange={(event) => {
            member.role = event.target.value
            setIsChanged(true)
            forceUpdate()
          }}
        >
          <MenuItem value={'Co-author'}>Co-author</MenuItem>
          <MenuItem value={'Reviewer'}>Reviewer</MenuItem>
        </Select>)
    }
  ]

  return <Datatable columns={columns} data={members}>
    <DatatableTable actions={DeleteAction} />
  </Datatable>
}

export const fetchUsers = (api, previousQuery, newQuery) => {
  return new Promise((resolve, reject) => {
    api.getUsers({prefix: newQuery})
      .then(users => {
        const withQueryInName = users.filter(user => user.name.toLowerCase().indexOf(newQuery) !== -1)
        withQueryInName.sort((a, b) => {
          const aValue = a.name.toLowerCase()
          const bValue = b.name.toLowerCase()
          if (aValue.startsWith(newQuery)) {
            return -1
          } else if (bValue.startsWith(newQuery)) {
            return 1
          } else {
            return 0
          }
        })
        resolve(withQueryInName.slice(0, 5))
      })
      .catch(err => {
        reject(err)
      })
  })
}

function AddMember({...props}) {
  const {api, raiseError} = props
  const [role, setRole] = useState('Co-author')
  const [suggestions, setSuggestions] = useState([])
  const [newMember, setNewMember] = useState([])
  const {members, setMembers, setIsChanged} = useContext(editMembersDialogContext)
  const [isDuplicated, setIsDuplicated] = useState(false)
  const [isValid, setIsValid] = useState(false)
  const [query, setQuery] = useState('')

  const handleInputChange = useCallback((event, value) => {
    const newQuery = value.toLowerCase()
    if (!(newQuery.startsWith(query) && suggestions.length === 0) || query === '') {
      fetchUsers(api, query, newQuery)
        .then(setSuggestions)
        .catch(err => {
          setSuggestions([])
          raiseError(err)
        })
    }
    setQuery(newQuery)
  }, [api, query, raiseError, suggestions.length])

  const debouncedHandleInputChange = useMemo(() => (
    debounce(handleInputChange, 700)
  ), [handleInputChange])

  const handleChange = useCallback((event, value) => {
    if (value && value?.user_id) {
      setNewMember(value)
      setIsValid(true)
      setIsDuplicated(members.map(member => member.user_id).includes(value.user_id))
    } else {
      setIsValid(false)
    }
  }, [members])

  const handleAdd = useCallback(() => {
    if (role) {
      if (!members.map(member => member.user_id).includes(newMember.user_id)) {
        newMember['role'] = role
        setMembers(members => [...members, newMember])
        setIsChanged(true)
      } else {
        setIsDuplicated(true)
      }
    }
  }, [members, newMember, role, setIsChanged, setMembers])

  return <React.Fragment>
    <AutoComplete
      style={{width: '100%'}}
      options={suggestions}
      getOptionLabel={option => (option.affiliation ? `${option.name} (${option.affiliation})` : option.name)}
      getOptionSelected={(option, value) => value ? option.user_id === value.user_id : false}
      onInputChange={debouncedHandleInputChange}
      onChange={handleChange}
      renderInput={params => (
        <TextField
          {...params}
          variant='filled'
          size='small'
          label='Search the name and select a user from the list'
          placeholder="Member's name"
          margin='normal'
          fullWidth
        />
      )}
    />
    <Box marginLeft={2}>
      <Typography hidden={!isDuplicated} color="error">
        The selected user is already in the members list
      </Typography>
    </Box>
    <FormControl variant='filled' size='small' fullWidth>
      <InputLabel htmlFor='role'>Select the member&apos;s role</InputLabel>
      <Select
        native
        onChange={(event) => setRole(event.target.value)}
        inputProps={{
          name: 'role',
          id: 'role'
        }}
      >
        <option value={'Co-author'}>Co-author</option>
        <option value={'Reviewer'}>Reviewer</option>
      </Select>
    </FormControl>
    <Box display="flex" justifyContent="end" paddingY={1}>
      <Button onClick={handleAdd} color="primary" variant="contained" disabled={isDuplicated || !isValid}>
        Add
      </Button>
    </Box>
  </React.Fragment>
}
AddMember.propTypes = {
  api: PropTypes.object.isRequired,
  raiseError: PropTypes.func.isRequired
}

const DeleteAction = React.memo((props) => {
  const {data} = props
  const {members, setMembers, setIsChanged} = useContext(editMembersDialogContext)

  const handleRemove = () => {
    const filteredMembers = members.filter(member => !(member.user_id === data.user_id))
    setMembers(filteredMembers)
    setIsChanged(true)
  }

  const isOwner = data.role === 'Main author'

  return <IconButton disabled={isOwner} onClick={handleRemove} data-testid='member-delete-button'>
    <Tooltip title="Remove the member">
      <DeleteIcon />
    </Tooltip>
  </IconButton>
})
DeleteAction.propTypes = {
  data: PropTypes.object.isRequired
}

function EditMembersDialog({...props}) {
  const classes = useStyles()
  const {disabled} = props
  const {api} = useApi()
  const {raiseError} = useErrors()
  const {uploadId, upload, updateUpload} = useUploadPageContext()
  const [open, setOpen] = useState(false)
  const [members, setMembers] = useState([])
  const [isChanged, setIsChanged] = useState(false)
  const [openConfirmDialog, setOpenConfirmDialog] = useState(false)
  const [isVisibleForAll, setIsVisibleForAll] = useState(false)

  const handleIsVisibleForAll = useCallback((event) => {
    setIsVisibleForAll(event.target.checked)
    setIsChanged(true)
  }, [])

  const getUsers = useCallback((user_ids, roles) => {
    return new Promise((resolve, reject) => {
      api.get(`users?user_id=${user_ids.join('&user_id=')}`)
        .then(response => {
          const members = response['data'].map((member, index) => {
            member.role = roles[index]
            return member
          })
          resolve(members)
        })
        .catch(error => reject(new Error('Unable to fetch the members' + error)))
      })
  }, [api])

  const fetchMembers = useCallback(() => {
    const user_ids = [upload.main_author].concat(upload.coauthors, upload.reviewers)
    const roles = ['Main author'].concat(upload.coauthors.map(_ => 'Co-author'), upload.reviewers.map(_ => 'Reviewer'))
    return getUsers(user_ids, roles)
  }, [getUsers, upload])

  const handleOpenDialog = () => {
    setMembers([])
    setIsChanged(false)
    setIsVisibleForAll(upload.reviewer_groups?.includes('all'))
    fetchMembers()
      .then(members => setMembers(members))
      .catch(error => raiseError(error))
    setOpen(true)
  }

  const handleDiscardChanges = () => {
    setOpenConfirmDialog(false)
    setOpen(false)
  }

  const handleSubmitChanges = () => {
    if (isChanged) {
      const newCoauthors = members.filter(member => member.role === 'Co-author').map(member => member.user_id)
      const newReviewers = members.filter(member => member.role === 'Reviewer').map(member => member.user_id)
      const allAction = isVisibleForAll ? 'add' : 'remove'
      api.post(`/uploads/${uploadId}/edit`, {
        'metadata': {
          'coauthors': newCoauthors,
          'reviewers': newReviewers,
          'reviewer_groups': {[allAction]: 'all'}
        }
      }).then(results => {
        updateUpload({upload: results.data})
        setOpen(false)
      }).catch(err => raiseError(err))
    } else {
      setOpen(false)
    }
  }

  const handleConfirm = () => {
    if (isChanged) {
      setOpenConfirmDialog(true)
    } else {
      setOpen(false)
    }
  }

  const contextValue = useMemo(() => ({
    members: members,
    setMembers: setMembers,
    isChanged: isChanged,
    setIsChanged: setIsChanged
  }), [members, setMembers, isChanged, setIsChanged])

  return <editMembersDialogContext.Provider value={contextValue}>
    <React.Fragment>
      <IconButton onClick={handleOpenDialog} disabled={disabled} data-testid='edit-members-action'>
        <Tooltip title="Manage upload members">
          <MembersIcon/>
        </Tooltip>
      </IconButton>
      {open && <Dialog classes={{paper: classes.dialog}} open={open} disableEscapeKeyDown data-testid='edit-members-dialog'>
        <DialogTitle>Manage upload members</DialogTitle>
        <DialogContent>
          <DialogContentText>
            You can add new members to this upload.
            <br/>
            The upload includes {upload?.entries} {upload?.entries === 1 ? 'entry' : 'entries'}.
          </DialogContentText>
          <Tooltip title="If checked, the upload can be viewed before publication even by unregistered users.">
            <FormControlLabel
              label="Visible for all"
              control={
                <Checkbox checked={isVisibleForAll} onChange={handleIsVisibleForAll} />
              }
            />
          </Tooltip>
          <Divider />
          <AddMember api={api} raiseError={raiseError} {...props}/>
          <MembersTable />
        </DialogContent>
        <DialogActions>
          <InviteUserDialog />
          <span style={{flexGrow: 1}} />
          <Button onClick={handleConfirm} color="secondary">
            Cancel
          </Button>
          <Button onClick={handleSubmitChanges} disabled={!isChanged} color="secondary">
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
  </editMembersDialogContext.Provider>
}
EditMembersDialog.propTypes = {
  disabled: PropTypes.bool
}

export default EditMembersDialog
