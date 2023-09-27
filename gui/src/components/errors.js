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
import React, { useContext } from 'react'
import PropTypes from 'prop-types'
import { SnackbarContent, IconButton, Snackbar, withStyles } from '@material-ui/core'
import CloseIcon from '@material-ui/icons/Close'

export class VersionMismatch extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'VersionMismatch'
  }
}

export const errorContext = React.createContext({
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
      if (error.name === 'BadRequest') {
        errorStr = error.message
      } else if (error.name === 'CannotReachApi') {
        errorStr = 'Cannot reach NOMAD, please try again later.'
      } else if (error.name === 'NotAuthorized') {
        errorStr = error.message
      } else if (error.name === 'DoesNotExist') {
        errorStr = 'You are trying to access information that does not exist. Please try again and let us know, if this error keeps happening.'
      } else if (error.name === 'VersionMismatch') {
        errorStr = 'There is a new GUI version available. Please press "shift" and reload the page.'
      } else if (error.message.startsWith('could not parse optimade')) {
        errorStr = 'The given OPTiMaDe query can not be parsed.'
      } else if (error.status === 422) {
        errorStr = `Unexpected api error: ${error.apiMessage[0].msg}. Please try again and let us know, if this error keeps happening.`
        console.log(error)
      } else if (error.message) {
        errorStr = `Unexpected error: "${error.message}". Please try again and let us know, if this error keeps happening.`
        console.log(error)
      }
    } else if (typeof error === 'string' || error instanceof String) {
      errorStr = `${error} Please close all NOMAD tabs and try again. Let us know, if this error keeps happening.`
    }

    if (this.state.errors.indexOf(errorStr) === -1) {
      this.setState({errors: [errorStr, ...this.state.errors]})
    }
  }

  onClose() {
    if (this.state.errors.length > 0) {
      this.setState({errors: this.state.errors.slice(1)})
    }
  }

  render() {
    const {children, classes} = this.props
    return (
      <errorContext.Provider value={this.state}>
        {children}
        <Snackbar
          anchorOrigin={{vertical: 'bottom', horizontal: 'left'}}
          open={this.state.errors.length > 0}
          onClose={this.handleClose}
        >
          <SnackbarContent
            className={classes.errorSnack}
            message={
              <span style={{color: 'white'}}>{'' + this.state.errors[0]}</span>
            }
            action={[
              <IconButton key={0} size="small" color="inherit" onClick={this.onClose.bind(this)}>
                <CloseIcon />
              </IconButton>
            ]}
          />
        </Snackbar>
      </errorContext.Provider>
    )
  }
}

export const ErrorSnacks = withStyles(ErrorSnacksUnstyled.styles)(ErrorSnacksUnstyled)

export function withErrors(Component) {
  function WithErrorComponent(props) {
    return (
      <errorContext.Consumer>
        {errorContext => <Component {...props} raiseError={errorContext.raiseError} />}
      </errorContext.Consumer>
    )
  }

  return WithErrorComponent
}

export class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props)
    this.state = { hasError: false }
  }

  static propTypes = {
    children: PropTypes.any,
    onError: PropTypes.func
  }

  static getDerivedStateFromError(_error) {
    // Update state so the next render will show the fallback UI.
    return { hasError: true }
  }

  componentDidCatch(error, errorInfo) {
    console.log('caught error in boundary', error, errorInfo)
    if (this.context) {
      this.context.raiseError('There has been a Javascript error.')
    }
  }

  render() {
    // If an error is caught, show the 'root' error message. This message will
    // be shown on top of everything does not contain any React or MUI specific
    // code to keep it functional no matter where the error originates from.
    if (this.state.hasError) {
      const error = document.getElementById("rootError")
      if (error) {
        error.style.display = 'block'
      }
      return null
    }

    return this.props.children
  }
}

ErrorBoundary.contextType = errorContext

export function useErrors() {
  return useContext(errorContext)
}
