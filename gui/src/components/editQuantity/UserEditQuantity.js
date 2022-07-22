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
import {useApi} from "../api"
import {useErrors} from "../errors"
import AutoComplete from "@material-ui/lab/Autocomplete"
import {debounce} from "lodash"
import {getFieldProps, TextFieldWithHelp} from './StringEditQuantity'
import {fetchUsers} from '../uploads/EditMembersDialog'

export const UserEditQuantity = React.memo((props) => {
  const {quantityDef, onChange, ...otherProps} = props
  const {helpDescription, ...otherFieldProps} = otherProps
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [query, setQuery] = useState('')
  const [user, setUser] = useState(undefined)
  const [suggestions, setSuggestions] = useState([])

  const handleInputChange = useCallback((event, value) => {
    if (!(event && event.type && event.type === 'change')) {
      return
    }
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
  }, [api, raiseError, query, suggestions])

  const getUser = useCallback((user_id) => {
    return new Promise((resolve, reject) => {
      if (!user_id) {
        resolve(undefined)
      }
      api.get(`users?user_id=${user_id}`)
        .then(response => {
          const user = response['data']?.[0]
          user ? resolve(user) : reject(new Error('Unable to find the member'))
        })
        .catch(error => {
          reject(new Error('Unable to find the member: ' + error))
        })
    })
  }, [api])

  useEffect(() => {
    getUser(otherProps.value).then(user => {
      setSuggestions(user ? [user] : [])
      setUser(user)
    })
  }, [getUser, otherProps.value])

  const debouncedHandleInputChange = useMemo(() => (
    debounce(handleInputChange, 700)
  ), [handleInputChange])

  const handleChange = useCallback((event, value) => {
    if (onChange) {
      setSuggestions(value ? [value] : [])
      setUser(value)
      onChange(value?.user_id === '' ? undefined : value?.user_id)
    }
  }, [onChange])

  return <React.Fragment>
    <AutoComplete
      style={{width: '100%'}}
      options={suggestions}
      getOptionLabel={option => option ? (option.affiliation ? `${option.name} (${option.affiliation})` : option.name) : ''}
      getOptionSelected={(option, value) => suggestions && value && option.user_id === value.user_id}
      onInputChange={debouncedHandleInputChange}
      onChange={handleChange}
      value={suggestions.includes(user) ? user : null}
      renderInput={params => (
        <TextFieldWithHelp
          {...params}
          variant='filled'
          size='small'
          placeholder="NOMAD member's name"
          margin='normal'
          fullWidth
          {...getFieldProps(quantityDef)}
          {...otherFieldProps}
        />
      )}
      data-testid='user-edit-quantity'
    />
  </React.Fragment>
})
UserEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}

export default UserEditQuantity
