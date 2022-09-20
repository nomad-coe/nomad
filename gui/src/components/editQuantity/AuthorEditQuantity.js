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
import React, {useCallback, useEffect, useState} from 'react'
import PropTypes from 'prop-types'
import {getFieldProps, TextFieldWithHelp} from './StringEditQuantity'
import {makeStyles, Typography} from '@material-ui/core'

const useStyles = makeStyles(theme => ({
  label: {
    marginLeft: theme.spacing(1)
  },
  fields: {
    marginTop: theme.spacing(1)
  }
}))

export const AuthorEditQuantity = React.memo((props) => {
  const {quantityDef, onChange, ...otherProps} = props
  const classes = useStyles()
  const [author, setAuthor] = useState(undefined)
  const [email, setEmail] = useState(otherProps?.value?.email)
  const [emailError, setEmailError] = useState(false)

  useEffect(() => {
    setAuthor(otherProps?.value)
  }, [otherProps?.value])

  const handleChange = useCallback((key, value) => {
    if (onChange) {
      const newValue = {...author}
      if (value) {
        newValue[key] = value
      } else {
        delete newValue[key]
      }
      onChange(newValue && Object.keys(newValue).length !== 0 ? newValue : undefined)
    }
  }, [author, onChange])

  const handleEmailChange = useCallback((value) => {
    setEmail(value)
    if (value.match(/^[^\s@]+@[^\s@]+\.[^\s@]+$/)) {
      handleChange('email', value)
      setEmailError(false)
    } else {
      setEmailError(true)
    }
  }, [handleChange])

  return <React.Fragment>
    <Typography className={classes.label} variant={'body1'}>
      {getFieldProps(quantityDef).label}
    </Typography>
    <TextFieldWithHelp
      className={classes.fields}
      label={'First name'}
      variant='filled'
      size='small'
      placeholder='First name'
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
      value={email}
      onChange={event => handleEmailChange(event.target.value)}
      error={emailError}
      helperText={emailError && 'The email is not valid!'}
      fullWidth
    />
    <TextFieldWithHelp
      className={classes.fields}
      label={'Address'}
      variant='filled'
      size='small'
      placeholder='Affiliation address'
      value={otherProps?.value?.affiliation_address}
      onChange={event => handleChange('affiliation_address', event.target.value)}
      fullWidth
    />
  </React.Fragment>
})
AuthorEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}

export default AuthorEditQuantity
