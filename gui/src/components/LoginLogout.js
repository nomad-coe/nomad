import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Typography from '@material-ui/core/Typography'
import { compose } from 'recompose'
import { Button, Link } from '@material-ui/core'
import { withApi } from './api'
import { keycloakBase, keycloakRealm } from '../config'

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
        marginRight: theme.spacing.unit * 2
      }
    },
    link: {
      color: 'white',
      textDecoration: 'underline'
    },
    button: {} // to allow overrides
  })

  render() {
    const { classes, variant, color, keycloak, user } = this.props

    if (keycloak.authenticated) {
      return (
        <div className={classes.root}>
          <Typography color="inherit" variant="body1">
            Welcome <Link
              className={classes.link}
              href={`${keycloakBase}/realms/${keycloakRealm}/account/`}>
              { user ? user.name : '...' }
            </Link>
          </Typography>
          <Button
            className={classes.button} style={{marginLeft: 8}}
            variant={variant} color={color}
            onClick={() => keycloak.logout()}
          >Logout</Button>
        </div>
      )
    } else {
      return (
        <div className={classes.root}>
          <Button
            className={classes.button} variant={variant} color={color} onClick={() => keycloak.login()}
          >Login / Register</Button>
        </div>
      )
    }
  }
}

export default compose(withApi(false), withStyles(LoginLogout.styles))(LoginLogout)
