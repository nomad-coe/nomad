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
import { SnackbarContent, IconButton, Snackbar, withStyles } from '@material-ui/core'
import ErrorIcon from '@material-ui/icons/Error'
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
    },
    icon: {
      fontSize: 20,
      opacity: 0.9,
      marginRight: theme.spacing(1)
    },
    message: {
      display: 'flex',
      alignItems: 'center'
    },
    errors: {
      marginLeft: theme.spacing(1)
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
        errorStr = 'Cannot reach NOMAD, please try again later.'
      } else if (error.name === 'NotAuthorized') {
        errorStr = error.message
      } else if (error.name === 'DoesNotExist') {
        errorStr = 'You are trying to access information that does not exist. Please try again and let us know, if this error keeps happening.'
      } else if (error.name === 'VersionMismatch') {
        errorStr = 'There is a new GUI version available. Please press "shift" and reload the page.'
      } else if (error.message.startsWith('could not parse optimade')) {
        errorStr = 'The given OPTiMaDe query can not be parsed.'
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
      <errorContext.Provider value={this.state}>
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
