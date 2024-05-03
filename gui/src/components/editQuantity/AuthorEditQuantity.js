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
import {getFieldProps, TextFieldWithHelp} from './StringEditQuantity'
import {Box, CircularProgress, FormLabel, InputAdornment, makeStyles, Typography} from '@material-ui/core'
import AutoComplete from '@material-ui/lab/Autocomplete'
import {debounce} from 'lodash'
import {fetchUsers} from '../uploads/EditMembersDialog'
import {useApi} from '../api'
import {useErrors} from '../errors'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

const useStyles = makeStyles(theme => ({
  label: {
    marginLeft: theme.spacing(1)
  },
  fields: {
    width: '100%',
    marginTop: theme.spacing(1)
  }
}))

export const AuthorEditQuantity = React.memo((props) => {
  const {quantityDef, onChange, ...otherProps} = props
  const {helpDescription, ...otherFieldProps} = otherProps
  const classes = useStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [author, setAuthor] = useState(undefined)
  const [user, setUser] = useState(undefined)
  const [email, setEmail] = useState(otherProps?.value?.email)
  const [query, setQuery] = useState('')
  const [suggestions, setSuggestions] = useState([])
  const [searching, setSearching] = useState(false)
  const userOnly = useMemo(() => quantityDef.type?.type_kind === 'User', [quantityDef])
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const searchUsers = useCallback((value) => {
    const newQuery = value.toLowerCase()
    if (!(newQuery.startsWith(query) && suggestions.length === 0) || query === '') {
      fetchUsers(api, query, newQuery)
        .then(setSuggestions)
        .catch(err => {
          raiseError(err)
          setSuggestions([])
        })
    }
    setQuery(newQuery)
    setSearching(false)
  }, [api, raiseError, query, suggestions])

  const handleChange = useCallback((key, value) => {
    if (onChange) {
      let newValue = {...author}
      if (key === 'user') {
        setSuggestions(value ? [value] : [])
        if (value) {
          newValue['user_id'] = value?.user_id
          newValue['first_name'] = value?.first_name
          newValue['last_name'] = value?.last_name
          newValue['affiliation'] = value?.affiliation
          newValue['affiliation_address'] = value?.affiliation_address
          delete newValue['email']
          setEmail(undefined)
          setAuthor(newValue)
        } else {
          newValue = undefined
        }
      } else {
        if (value) {
          newValue[key] = value
        } else {
          delete newValue[key]
        }
      }
      onChange(newValue && Object.keys(newValue).length !== 0 ? newValue : undefined)
    }
  }, [author, onChange])

  useEffect(() => {
    const initUser = async (user_id) => {
      if (user_id) {
        const users = await api.getUsers({user_id: user_id})
        if (users.length > 0) {
          setSuggestions(users)
          setUser(users[0])
        } else {
          setSuggestions([])
          setUser(undefined)
        }
      } else {
        setUser(undefined)
        setSuggestions([])
      }
    }

    initUser(otherFieldProps.value?.user_id)

    if (typeof otherProps?.value === 'string' && otherProps.value.includes('@')) {
      api.getUsers({email: otherProps.value})
        .then(users => {
          if (users.length > 0) {
            handleChange('user', users[0])
            setSuggestions([users[0]])
          } else {
            setSuggestions([])
          }
        })
    } else {
      setAuthor(otherProps?.value)
    }
  }, [api, handleChange, otherFieldProps.value, otherProps.value, quantityDef.default])

  const isValidEmail = useCallback((value) => value ? value.match(/^[^\s@]+@[^\s@]+\.[^\s@]+$/) : true, [])

  const handleEmailChange = useCallback((value) => {
    setEmail(value)
    if (isValidEmail(value)) {
      handleChange('email', value)
    }
  }, [handleChange, isValidEmail])

  const debouncedSearchUsers = useMemo(() => debounce(searchUsers, 700), [searchUsers])

  const handleUserInputChange = useCallback((event, value) => {
    if (!(event && event.type && event.type === 'change')) {
      return
    }
    setSearching(true)
    debouncedSearchUsers(value)
  }, [debouncedSearchUsers])

  const userField = (
    <AutoComplete
      options={suggestions}
      getOptionLabel={option => option ? (option.affiliation ? `${option.name} (${option.affiliation})` : option.name) : ''}
      getOptionSelected={(option, value) => suggestions && value && value.user_id && option.user_id === value.user_id}
      onInputChange={handleUserInputChange}
      onChange={(event, value) => handleChange('user', value)}
      value={user && suggestions.some(option => option.user_id === user.user_id) ? user : null}
      renderInput={params => (
        <TextFieldWithHelp
          className={classes.fields}
          {...params}
          InputProps={(searching ? {
            endAdornment: (searching && <InputAdornment position="end">
              <CircularProgress color="inherit" size={20} />
            </InputAdornment>)
          } : params.InputProps)}
          variant='filled'
          size='small'
          placeholder="Search member's name"
          margin='normal'
          fullWidth
          {...getFieldProps(quantityDef)}
          label={(userOnly ? label : 'User account')}
        />
      )}
      data-testid='user-edit-quantity'
    />
  )

  return <React.Fragment>
    {userOnly ? userField : (
      <Box marginTop={2} marginBottom={2}>
        <FormLabel component="legend">
          {label}
        </FormLabel>
        {userField}
        <Typography className={classes.label} variant={'caption'}>
          {'Or provide author information'}
        </Typography>
        <TextFieldWithHelp
          label={'First name'}
          variant='filled'
          size='small'
          placeholder='First name'
          inputProps={{readOnly: !!otherProps?.value?.user_id}}
          value={otherProps?.value?.first_name}
          onChange={event => handleChange('first_name', event.target.value)}
          fullWidth
        />
        <TextFieldWithHelp
          className={classes.fields}
          label={'Last name'}
          variant='filled'
          size='small'
          placeholder='Last name'
          inputProps={{readOnly: !!otherProps?.value?.user_id}}
          value={otherProps?.value?.last_name}
          onChange={event => handleChange('last_name', event.target.value)}
          fullWidth
        />
        <TextFieldWithHelp
          className={classes.fields}
          label={'Affiliation'}
          variant='filled'
          size='small'
          placeholder='Affiliation'
          inputProps={{readOnly: !!otherProps?.value?.user_id}}
          value={otherProps?.value?.affiliation}
          onChange={event => handleChange('affiliation', event.target.value)}
          fullWidth
        />
        <TextFieldWithHelp
          className={classes.fields}
          label={'Email'}
          variant='filled'
          size='small'
          placeholder="Email address"
          inputProps={{readOnly: !!otherProps?.value?.user_id}}
          value={email}
          onChange={event => handleEmailChange(event.target.value)}
          error={!isValidEmail(email)}
          helperText={!isValidEmail(email) && 'The email is not valid!'}
          fullWidth
        />
        <TextFieldWithHelp
          className={classes.fields}
          label={'Address'}
          variant='filled'
          size='small'
          placeholder='Affiliation address'
          inputProps={{readOnly: !!otherProps?.value?.user_id}}
          value={otherProps?.value?.affiliation_address}
          onChange={event => handleChange('affiliation_address', event.target.value)}
          fullWidth
        />
      </Box>
    )}
  </React.Fragment>
})
AuthorEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.object,
  onChange: PropTypes.func
}

export default AuthorEditQuantity
