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
import React from 'react'
import PropTypes from 'prop-types'
import Typography from '@material-ui/core/Typography'
import { Button, Link, makeStyles } from '@material-ui/core'
import { keycloakBase, keycloakRealm } from '../config'
import LoginIcon from '@material-ui/icons/AccountCircle'
import { useApi } from './api'
import { useKeycloak } from '@react-keycloak/web'

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    '& p': {
      marginRight: theme.spacing(2)
    }
  },
  link: {
    textDecoration: 'underline'
  },
  button: {} // to allow overrides
}))

const LoginLogout = React.memo(function LoginLogout(props) {
  const classes = useStyles()
  const {user} = useApi()
  const {keycloak} = useKeycloak()
  const {variant, color} = props

  if (keycloak?.authenticated) {
    return (
      <div className={classes.root}>
        <Typography color="primary" variant="body1">
          Welcome <Link
            className={classes.link}
            href={`${keycloakBase.replace(/\/$/, '')}/realms/${keycloakRealm}/account/`}>
            { user ? user.name : '...' }
          </Link>
        </Typography>
        <Button
          className={classes.button} style={{marginLeft: 8}}
          variant={variant} color={color}
          onClick={() => keycloak.logout()}
          startIcon={<LoginIcon/>}
        >Logout</Button>
      </div>
    )
  } else {
    return (
      <div className={classes.root}>
        <Button
          className={classes.button} variant={variant} color={color}
          startIcon={<LoginIcon/>}
          onClick={() => keycloak.login()}
          data-testid='login-register-button'
        >Login / Register</Button>
      </div>
    )
  }
})
LoginLogout.propTypes = {
  variant: PropTypes.string,
  color: PropTypes.string
}

export default LoginLogout
