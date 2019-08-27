import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Typography from '@material-ui/core/Typography'
import { compose } from 'recompose'
import { Button } from '@material-ui/core'
import { withApi } from './api'

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
    button: {} // to allow overrides
  })

  render() {
    const { classes, variant, color, keycloak, user } = this.props

    if (keycloak.authenticated) {
      return (
        <div className={classes.root}>
          <Typography color="inherit" variant="body1">
            Welcome { user ? user.name : '...' }
          </Typography>
          <Button
            className={classes.button}
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
          >Login</Button>
        </div>
      )
    }
  }
}

export default compose(withApi(false), withStyles(LoginLogout.styles))(LoginLogout)
