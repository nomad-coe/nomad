import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Typography from '@material-ui/core/Typography'
import { compose } from 'recompose'
import { Button, Link } from '@material-ui/core'
import { withApi } from './api'
import { keycloakBase, keycloakRealm } from '../config'
import LoginIcon from '@material-ui/icons/AccountCircle'

class LoginLogout extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    variant: PropTypes.string,
    color: PropTypes.string,
    user: PropTypes.object,
    keycloak: PropTypes.object.isRequired
  }

  static styles = theme => ({
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
  })

  render() {
    const { classes, variant, color, keycloak, user } = this.props

    if (keycloak.authenticated) {
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
          >Login / Register</Button>
        </div>
      )
    }
  }
}

export default compose(withApi(false), withStyles(LoginLogout.styles))(LoginLogout)
