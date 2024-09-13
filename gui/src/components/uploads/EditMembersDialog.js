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
  FormGroup,
  Dialog,
  DialogContent,
  DialogTitle,
  IconButton,
  makeStyles,
  MenuItem,
  Select,
  TextField,
  Tooltip
} from '@material-ui/core'
import Button from '@material-ui/core/Button'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContentText from '@material-ui/core/DialogContentText'
import DeleteIcon from '@material-ui/icons/Delete'
import AutoComplete from '@material-ui/lab/Autocomplete'
import _, { debounce } from 'lodash'
import PropTypes from 'prop-types'
import React, { useCallback, useContext, useEffect, useMemo, useReducer, useState } from 'react'
import { isUploadVisibleForAll } from '../../utils'
import { useApi } from '../api'
import { Datatable, DatatableTable } from '../datatable/Datatable'
import { useErrors } from '../errors'
import {
  fetchGroupsByIds,
  fetchGroupsByNameSearch,
  fetchUsersByIds,
  fetchUsersByNamePrefix
} from '../utils/apiUtils'
import { InviteUserDialog } from './InviteUserDialog'
import { useUploadPageContext } from './UploadPageContext'
import {
  getMemberIdsByTypeAndRole,
  groupsToMembers,
  groupToMember,
  MemberIcon,
  ROLES,
  TYPES,
  userToMember
} from './memberUtils'
import { uploadMembersGroupSearchEnabled } from '../../config'

export const editMembersDialogContext = React.createContext()

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%'
  },
  memberSearch: {
    columnGap: theme.spacing(1)
  },
  memberSearchType: {
    width: '7em'
  },
  memberSearchText: {
    flexGrow: 1
  }
}))

function MembersTable() {
  const {members, setIsChanged} = useContext(editMembersDialogContext)
  const forceUpdate = useReducer(bool => !bool)[1]

  const columns = [
    {key: '', align: 'right', render: member => (<MemberIcon type={member.type} />)},
    {key: 'Name', align: 'left', render: member => member.name},
    {key: 'Affiliation', align: 'left', render: member => member.affiliation},
    {
      key: 'Role',
      align: 'left',
      render: member => (member.role === ROLES.MAIN_AUTHOR ? member.role
        : <Select
          value={member.role}
          onChange={(event) => {
            member.role = event.target.value
            setIsChanged(true)
            forceUpdate()
          }}
        >
          <MenuItem value={ROLES.COAUTHOR}>Co-author</MenuItem>
          <MenuItem value={ROLES.REVIEWER}>Reviewer</MenuItem>
        </Select>)
    }
  ]

  return <Datatable columns={columns} data={members}>
    <DatatableTable actions={DeleteAction} />
  </Datatable>
}

function AddMember() {
  const classes = useStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [suggestions, setSuggestions] = useState([])
  const {members, setMembers, setIsChanged} = useContext(editMembersDialogContext)
  const [query, setQuery] = useState('')
  const [searchType, setSearchType] = useState(TYPES.USER)
  const [selectedSuggestion, setSelectedSuggestion] = useState(null)
  const [suggestionsAreLoading, setSuggestionsAreLoading] = useState(false)

  const handleInputChange = useCallback(async (event, value) => {
    const newQuery = value.toLowerCase()
    const oldQuery = query

    if (oldQuery !== '' && newQuery.startsWith(oldQuery) && suggestions.length === 0) {
      setQuery(newQuery)
      return
    }

    try {
      setSuggestionsAreLoading(true)
      let newSuggestions = []
      if (searchType === TYPES.GROUP) {
        const groups = await fetchGroupsByNameSearch(api, newQuery)
        newSuggestions = await groupsToMembers(groups, api)
      } else {
        const users = await fetchUsersByNamePrefix(api, newQuery)
        newSuggestions = users.map(user => userToMember(user))
      }

      const filteredSuggestions = newSuggestions.filter(
        suggestion => !members.some(member => member.id === suggestion.id)
      )

      setSuggestions(filteredSuggestions)
    } catch (err) {
      setSuggestions([])
      raiseError(err)
    } finally {
      setSuggestionsAreLoading(false)
    }

    setQuery(newQuery)
  }, [api, searchType, members, query, raiseError, suggestions.length])

  const debouncedHandleInputChange = useMemo(() => (
    debounce(handleInputChange, 700)
  ), [handleInputChange])

  const handleSearchTypeChange = useCallback((event) => {
    setSearchType(event.target.value)
    setSelectedSuggestion(null)
  }, [])

  const handleSearchSelectionChange = useCallback((event, value) => {
    setSelectedSuggestion(null)

    if (members.some(member => member.id === value?.id)) {
      return
    }

    setMembers([...members, value])
    setIsChanged(true)
  }, [members, setMembers, setIsChanged])

  return <React.Fragment>
    <FormGroup row className={classes.memberSearch}>
      {uploadMembersGroupSearchEnabled && <TextField
        select
        className={classes.memberSearchType}
        id='search-type'
        size='small'
        margin='normal'
        variant='filled'
        value={searchType}
        onChange={handleSearchTypeChange}
      >
        <MenuItem value={TYPES.USER}>User</MenuItem>
        <MenuItem value={TYPES.GROUP}>Group</MenuItem>
      </TextField>}
      <AutoComplete
        className={classes.memberSearchText}
        options={suggestions}
        getOptionLabel={option => (option.affiliation ? `${option.name} (${option.affiliation})` : option.name)}
        getOptionSelected={(option, value) => value ? option.id === value.id : false}
        onInputChange={debouncedHandleInputChange}
        onChange={handleSearchSelectionChange}
        value={selectedSuggestion}
        blurOnSelect={true}
        loading={suggestionsAreLoading}
        renderInput={params => (
          <TextField
            {...params}
            size='small'
            margin='normal'
            variant='filled'
            label={`Search the name and select a ${searchType} from the list`}
            placeholder={searchType === TYPES.GROUP ? "Group's name" : "User's name"}
          />
        )}
      />
    </FormGroup>
  </React.Fragment>
}

const DeleteAction = React.memo((props) => {
  const {data} = props
  const {members, setMembers, setIsChanged} = useContext(editMembersDialogContext)

  const handleRemove = () => {
    const filteredMembers = members.filter(member => member.id !== data.id)
    setMembers(filteredMembers)
    setIsChanged(true)
  }

  const isOwner = data.role === ROLES.MAIN_AUTHOR

  return <IconButton disabled={isOwner} onClick={handleRemove} data-testid='member-delete-button'>
    <Tooltip title="Remove the member">
      <DeleteIcon />
    </Tooltip>
  </IconButton>
})
DeleteAction.propTypes = {
  data: PropTypes.object.isRequired
}

const EditMembersDialog = ({open, setOpen}) => {
  const classes = useStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const {uploadId, upload, updateUpload} = useUploadPageContext()
  const [members, setMembers] = useState([])
  const [isChanged, setIsChanged] = useState(false)
  const [openConfirmDialog, setOpenConfirmDialog] = useState(false)

  const fetchMembers = useCallback(async () => {
    const { main_author, coauthors, reviewers, coauthor_groups } = upload
    const reviewer_groups = upload.reviewer_groups.filter(group => group !== 'all')

    const groupIds = [...coauthor_groups, ...reviewer_groups]
    const groupRoles = [...coauthor_groups.map(_ => ROLES.COAUTHOR),
                        ...reviewer_groups.map(_ => ROLES.REVIEWER)]
    const groups = await fetchGroupsByIds(api, groupIds)
    const ownerIds = groups.map(group => group.owner)
    const groupsDict = _.keyBy(groups, 'group_id')

    const memberUserIds = [main_author, ...coauthors, ...reviewers]
    const userRoles = [ROLES.MAIN_AUTHOR,
                       ...coauthors.map(_ => ROLES.COAUTHOR),
                       ...reviewers.map(_ => ROLES.REVIEWER)]
    const users = await fetchUsersByIds(api, [...memberUserIds, ...ownerIds])
    const usersDict = _.keyBy(users, 'user_id')

    const memberUsers = memberUserIds.map((userId, index) =>
      userToMember(usersDict[userId], userRoles[index])
    )
    const memberGroups = groupIds.map((groupId, index) =>
      groupToMember(groupsDict[groupId], usersDict, groupRoles[index])
    )

    return [...memberUsers, ...memberGroups]
  }, [api, upload])

  useEffect(() => {
    if (!open) {
      return
    }

    setMembers([])
    setIsChanged(false)
    fetchMembers()
      .then(setMembers)
      .catch(raiseError)
  }, [fetchMembers, open, raiseError])

  const handleDiscardChanges = () => {
    setOpenConfirmDialog(false)
    setOpen(false)
  }

  const handleSubmitChanges = () => {
    if (!isChanged) {
      setOpen(false)
      return
    }

    const metadata = {
      'coauthors': getMemberIdsByTypeAndRole(members, TYPES.USER, ROLES.COAUTHOR),
      'reviewers': getMemberIdsByTypeAndRole(members, TYPES.USER, ROLES.REVIEWER),
      'coauthor_groups': getMemberIdsByTypeAndRole(members, TYPES.GROUP, ROLES.COAUTHOR),
      'reviewer_groups': getMemberIdsByTypeAndRole(members, TYPES.GROUP, ROLES.REVIEWER)
    }

    if (isUploadVisibleForAll(upload)) {
      metadata.reviewer_groups.push('all')
    }

    api.post(`/uploads/${uploadId}/edit`, {'metadata': metadata})
      .then(results => {
        updateUpload({upload: results.data})
        setOpen(false)
      })
      .catch(raiseError)
  }

  const handleCancel = () => {
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
      {open && <Dialog classes={{paper: classes.dialog}} open={open} disableEscapeKeyDown data-testid='edit-members-dialog'>
        <DialogTitle>Edit upload members</DialogTitle>
        <DialogContent>
          <DialogContentText>
            You can add and remove members for this upload and change their access.
            Reviewers can view the upload, co-authors can also edit it.
          </DialogContentText>
          <AddMember />
          <MembersTable />
        </DialogContent>
        <DialogActions>
          <InviteUserDialog />
          <span style={{flexGrow: 1}} />
          <Button onClick={handleCancel} color="secondary">
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
  open: PropTypes.bool.isRequired,
  setOpen: PropTypes.func.isRequired
}

export default EditMembersDialog
