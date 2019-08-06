import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Typography from '@material-ui/core/Typography'
import { compose } from 'recompose'
import { Button, DialogTitle, DialogContent, DialogContentText, TextField, DialogActions,
  Dialog, FormGroup } from '@material-ui/core'
import { withApi } from './api'

class LoginLogout extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    isLoggingIn: PropTypes.bool,
    user: PropTypes.object,
    login: PropTypes.func.isRequired,
    logout: PropTypes.func.isRequired,
    variant: PropTypes.string,
    color: PropTypes.string,
    onLoggedIn: PropTypes.func,
    onLoggedOut: PropTypes.func
  }

  static styles = theme => ({
    root: {
      display: 'flex',
      alignItems: 'center',
      '& p': {
        marginRight: theme.spacing.unit * 2
      }
    },
    button: {}, // to allow overrides
    buttonDisabled: {},
    errorText: {
      marginTop: theme.spacing.unit,
      marginBottom: theme.spacing.unit
    }
  })

  constructor(props) {
    super(props)
    this.handleLogout = this.handleLogout.bind(this)
    this.handleChange = this.handleChange.bind(this)
    this.handleKeyPress = this.handleKeyPress.bind(this)
  }

  state = {
    loginDialogOpen: false,
    userName: '',
    password: '',
    failure: false
  }

  componentDidMount() {
    this._ismounted = true
  }

  componentWillUnmount() {
    this._ismounted = false
  }

  handleLoginDialogClosed(withLogin) {
    if (withLogin) {
      this.props.login(this.state.userName, this.state.password, (success) => {
        if (success && this.props.onLoggedIn) {
          this.props.onLoggedIn()
        }

        if (this._ismounted) {
          if (success) {
            this.setState({loginDialogOpen: false, failure: false})
          } else {
            this.setState({failure: true, loginDialogOpen: true})
          }
        }
      })
    } else {
      if (this._ismounted) {
        this.setState({failure: false, userName: '', password: '', loginDialogOpen: false})
      }
    }
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value, failure: false
    })
  }

  handleLogout() {
    this.props.logout()
    if (this.props.onLoggedOut) {
      this.props.onLoggedOut()
    }
  }

  handleKeyPress(ev) {
    if (ev.key === 'Enter') {
      ev.preventDefault()
      this.handleLoginDialogClosed(true)
    }
  }

  render() {
    const { classes, user, variant, color, isLoggingIn } = this.props
    const { failure } = this.state
    if (user) {
      return (
        <div className={classes.root}>
          <Typography color="inherit" variant="body1">
            Welcome, {user.first_name} {user.last_name}
          </Typography>
          <Button
            className={classes.button}
            variant={variant} color={color}
            onClick={this.handleLogout}
          >Logout</Button>
        </div>
      )
    } else {
      return (
        <div className={classes.root}>
          <Button
            className={isLoggingIn ? classes.buttonDisabled : classes.button} variant={variant} color={color} disabled={isLoggingIn}
            onClick={() => this.setState({loginDialogOpen: true})}
          >Login</Button>
          <Dialog
            disableBackdropClick disableEscapeKeyDown
            open={this.state.loginDialogOpen}
            onClose={() => this.handleLoginDialogClosed(false)}
          >
            <DialogTitle>Login</DialogTitle>
            <DialogContent>
              <DialogContentText>
                To login, please enter your email address and password. If you
                do not have an account, please go to the nomad repository and
                create one.
              </DialogContentText>
              {failure ? <DialogContentText className={classes.errorText} color="error">Wrong username or password!</DialogContentText> : ''}
              <form>
                <FormGroup>
                  <TextField
                    autoComplete="username"
                    disabled={isLoggingIn}
                    autoFocus
                    margin="dense"
                    id="uaseName"
                    label="Email Address"
                    type="email"
                    fullWidth
                    value={this.state.userName}
                    onChange={this.handleChange('userName')}
                    onKeyPress={this.handleKeyPress}
                  />
                  <TextField
                    autoComplete="current-password"
                    disabled={isLoggingIn}
                    margin="dense"
                    id="password"
                    label="Password"
                    type="password"
                    fullWidth
                    value={this.state.password}
                    onChange={this.handleChange('password')}
                    onKeyPress={this.handleKeyPress}
                  />
                </FormGroup>
              </form>
            </DialogContent>
            <DialogActions>
              <Button onClick={() => this.handleLoginDialogClosed(false)} color="primary">
                Cancel
              </Button>
              <Button onClick={() => this.handleLoginDialogClosed(true)} color="primary"
                disabled={this.state.userName === '' || this.state.password === ''}
              >
                Login
              </Button>
            </DialogActions>
          </Dialog>
        </div>
      )
    }
  }
}

export default compose(withApi(false), withStyles(LoginLogout.styles))(LoginLogout)
