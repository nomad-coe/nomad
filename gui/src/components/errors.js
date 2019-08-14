import React from 'react'
import PropTypes from 'prop-types'
import { SnackbarContent, IconButton, Snackbar, withStyles } from '@material-ui/core'
import ErrorIcon from '@material-ui/icons/Error'
import CloseIcon from '@material-ui/icons/Close'

export const ErrorContext = React.createContext({
  errors: [],
  raiseError: () => { throw Error('Error context used incorrectly.') }
})

class ErrorSnacksUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any.isRequired
  }
  static styles = theme => ({
    root: {},
    errorSnack: {
      backgroundColor: theme.palette.error.dark
    },
    icon: {
      fontSize: 20,
      opacity: 0.9,
      marginRight: theme.spacing.unit
    },
    message: {
      display: 'flex',
      alignItems: 'center'
    },
    errors: {
      marginLeft: theme.spacing.unit
    }
  })

  state = {
    errors: [],
    raiseError: (error) => {
      this.onError(error)
    }
  }

  onError(error) {
    let errorStr = 'Unexpected error. Please try again and let us know, if this error keeps happening.'
    if (error instanceof Error) {
      if (error.name === 'CannotReachApi') {
        errorStr = 'Cannot reach the NOMAD API, please try again later.'
      } else if (error.name === 'DoesNotExist') {
        errorStr = 'You are trying to access information that does not exist. Please try again and let us know, if this error keeps happening.'
      } else if (error.message) {
        errorStr = `Unexpected error: "${error.message}". Please try again and let us know, if this error keeps happening.`
      }
    } else if (error instanceof String) {
      errorStr = `Unexpected error: "${error}". Please try again and let us know, if this error keeps happening.`
    }

    if (this.state.errors.indexOf(errorStr) === -1) {
      this.setState({errors: [errorStr, ...this.state.errors]})
    }
  }

  onClose() {
    this.setState({errors: []})
  }

  render() {
    const {children, classes} = this.props
    return (
      <ErrorContext.Provider value={this.state}>
        {children}
        <Snackbar
          anchorOrigin={{vertical: 'bottom', horizontal: 'left'}}
          open={this.state.errors.length > 0}
          // autoHideDuration={6000}
          onClose={this.handleClose}
        >
          <SnackbarContent
            className={classes.errorSnack}
            message={
              <span id="client-snackbar" className={classes.message}>
                <ErrorIcon className={classes.icon} />
                <div className={classes.errors}>
                  {this.state.errors.map((error, index) => (
                    <p key={index}>{'' + error}</p>
                  ))}
                </div>
              </span>
            }
            action={[
              <IconButton key={0} color="inherit" onClick={this.onClose.bind(this)}>
                <CloseIcon />
              </IconButton>
            ]}
          />
        </Snackbar>
      </ErrorContext.Provider>
    )
  }
}

export const ErrorSnacks = withStyles(ErrorSnacksUnstyled.styles)(ErrorSnacksUnstyled)

export function withErrors(Component) {
  function WithErrorComponent(props) {
    return (
      <ErrorContext.Consumer>
        {errorContext => <Component {...props} raiseError={errorContext.raiseError} />}
      </ErrorContext.Consumer>
    )
  }

  return WithErrorComponent
}
