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
import React, { useState } from 'react'
import PropTypes from 'prop-types'
import Alert from '@material-ui/lab/Alert'
import { hasWebGLSupport } from '../utils'

export class ErrorHandler extends React.Component {
  state = {
    hasError: false
  }

  static getDerivedStateFromError(error) {
    return {
      hasError: true,
      error: error
    }
  }

  componentDidCatch(error, errorInfo) {
    console.log(error, errorInfo)
  }

  render() {
    if (this.state.hasError) {
      const msg = typeof this.props.message === 'string' ? this.props.message : this.props.message(this.state.error)
      return <Alert
        severity="error"
        className={this.props.className}
        classes={this.props.classes}
      >
        {msg}
      </Alert>
    }
    return this.props.children
  }
}
ErrorHandler.propTypes = ({
  children: PropTypes.node,
  message: PropTypes.oneOfType([PropTypes.string, PropTypes.func]), // Provide either a fixed error message or a callback that will receive the error details.
  classes: PropTypes.object,
  className: PropTypes.string
})

export const withErrorHandler = (message) => (WrappedComponent) => props => {
  return <ErrorHandler message={message}>
    <WrappedComponent {...props}></WrappedComponent>
  </ErrorHandler>
}
withErrorHandler.propTypes = ({
  message: PropTypes.oneOfType([PropTypes.string, PropTypes.func]) // Provide either a fixed error message or a callback that will receive the error details.
})

export const webGlError = 'Could not display the visualization as your browser does not support WebGL content.'
export const withWebGLErrorHandler = WrappedComponent => props => {
  const hasWebGL = useState(hasWebGLSupport())[0]

  // If WebGL is not available, the content cannot be shown.
  if (hasWebGL) {
    return WrappedComponent({...props})
  } else {
    return <Alert
      severity="info"
    >
      {webGlError}
    </Alert>
  }
}

withErrorHandler.propTypes = ({
  message: PropTypes.oneOfType([PropTypes.string, PropTypes.func]) // Provide either a fixed error message or a callback that will receive the error details.
})
