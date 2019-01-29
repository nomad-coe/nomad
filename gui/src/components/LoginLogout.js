import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Typography from '@material-ui/core/Typography'
import { compose } from 'recompose'
import { Button, DialogTitle, DialogContent, DialogContentText, TextField, DialogActions,
  Dialog, FormGroup, LinearProgress } from '@material-ui/core'
import { withApi } from './api'

class LoginLogout extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    login: PropTypes.func.isRequired,
    logout: PropTypes.func.isRequired,
    variant: PropTypes.string,
    color: PropTypes.string
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
    errorText: {
      marginTop: theme.spacing.unit,
      marginBottom: theme.spacing.unit
    }
  })

  constructor(props) {
    super(props)
    this.handleLogout = this.handleLogout.bind(this)
    this.handleChange = this.handleChange.bind(this)
  }

  state = {
    loginDialogOpen: false,
    userName: '',
    password: '',
    loggingIn: false,
    failure: false
  }

  componentDidMount() {
    this._ismounted = true
  }

  componentWillUnmount() {
    this._ismounted = false
  }

  handleLoginDialogClosed(withLogin) {
    this.setState({loginDialogOpen: false})
    if (withLogin) {
      if (this._ismounted) {
        this.setState({loggingIn: true})
      }
      this.props.login(this.state.userName, this.state.password, (success) => {
        if (this._ismounted) {
          if (success) {
            this.setState({loggingIn: false, loginDialogOpen: false, failure: false})
          } else {
            this.setState({loggingIn: false, failure: true, loginDialogOpen: true})
          }
        }
      })
    } else {
      if (this._ismounted) {
        this.setState({loggingIn: false, failure: false, userName: '', password: '', loginDialogOpen: false})
      }
    }
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    })
  }

  handleLogout() {
    this.props.logout()
  }

  render() {
    const { classes, user, variant, color } = this.props
    const { loggingIn, failure } = this.state
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
            className={classes.button} variant={variant} color={color}
            onClick={() => this.setState({loginDialogOpen: true})}
          >Login</Button>
          <Dialog
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
              {loggingIn ? <LinearProgress/> : ''}
              {failure ? <DialogContentText className={classes.errorText} color="error">Wrong username or password!</DialogContentText> : ''}
              <FormGroup>
                <TextField
                  disabled={loggingIn}
                  autoFocus
                  margin="dense"
                  id="uaseName"
                  label="Email Address"
                  type="email"
                  fullWidth
                  value={this.state.userName}
                  onChange={this.handleChange('userName')}
                />
                <TextField
                  disabled={loggingIn}
                  margin="dense"
                  id="password"
                  label="Password"
                  type="password"
                  fullWidth
                  value={this.state.password}
                  onChange={this.handleChange('password')}
                />
              </FormGroup>
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
