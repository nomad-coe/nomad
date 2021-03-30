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
import React, { useContext, useEffect, useCallback, useRef } from 'react'
import PropTypes from 'prop-types'
import { withErrors } from './errors'
import { apiBase } from '../config'
import { Typography, withStyles } from '@material-ui/core'
import LoginLogout from './LoginLogout'
import { compose } from 'recompose'
import { withKeycloak } from 'react-keycloak'
import axios from 'axios'

export const apiContext = React.createContext()

export class DoesNotExist extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'DoesNotExist'
  }
}

export class NotAuthorized extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'NotAuthorized'
  }
}

export class ApiError extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'CannotReachApi'
  }
}

function handleApiError(e) {
  if (e.name === 'CannotReachApi' || e.name === 'NotAuthorized' || e.name === 'DoesNotExist') {
    throw e
  }

  let error = null
  if (e.response) {
    const body = e.response.body
    const message = (body && (body.message || body.description)) || e.response.statusText
    const errorMessage = `${message} (${e.response.status})`
    if (e.response.status === 404) {
      error = new DoesNotExist(errorMessage)
    } else if (e.response.status === 401) {
      error = new NotAuthorized(errorMessage)
    } else if (e.response.status === 502) {
      error = new ApiError(errorMessage)
    } else {
      error = new Error(errorMessage)
    }
    error.status = e.response.status
    error.apiMessage = message
  } else {
    if (e.message === 'Failed to fetch') {
      error = new ApiError(e.message)
      error.status = 400
    } else {
      const errorMessage = e.status ? `${e} (${e.status})` : '' + e
      error = new Error(errorMessage)
    }
  }
  throw error
}

class Api {
  constructor(keycloak) {
    this.keycloak = keycloak
    this.axios = axios.create({
      baseURL: `${apiBase}/v1`
    })

    this.loadingHandler = []
    this.loading = 0

    this.onFinishLoading = () => {
      this.loading--
      this.loadingHandler.forEach(handler => handler(this.loading))
    }
    this.onStartLoading = () => {
      this.loading++
      this.loadingHandler.forEach(handler => handler(this.loading))
    }

    Api.uploadIds = 0
  }

  /**
   * Fetches the up-to-date authorization headers. Performs a token update if
   * the current token has expired.
   * @returns Object containing the authorization header.
   */
  async authHeaders() {
    if (this.keycloak.token) {
      return new Promise((resolve, reject) => {
        this.keycloak.updateToken()
          .success(() => {
            resolve({
              'Authorization': `Bearer ${this.keycloak.token}`
            })
          })
          .error(() => {
            reject(new ApiError())
          })
      })
    } else {
      return {}
    }
  }

  /**
   * Returns the entry information that is stored in the search index.
   * @param {string} entryId
   * @returns Object containing the search index contents.
   */
  async entry(entryId) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const entry = await this.axios.get(`/entries/${entryId}`, auth)
      return entry.data
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }

  /**
   * Returns section_results from the archive corresponding to the given entry.
   * All references within the section are resolved by the server before
   * sending.
   * @param {string} entryId
   * @returns Object containing section_results
   */
  async results(entryId) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const entry = await this.axios.post(
        `/entries/archive/query`,
        {
          query: {calc_id: entryId},
          pagination: {size: 1},
          required: {results: '*'}
        },
        auth
      )
      return parse(entry.data)
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }

  /**
   * Return the raw file metadata for a given entry.
   * @param {string} entryId
   * @returns Object containing the raw file metadata.
   */
  async getRawFileListFromCalc(entryId) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const entry = await this.axios.get(`/entries/${entryId}/raw`, auth)
      return entry.data
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }
}

/**
 * Custom JSON parse function that can handle NaNs that can be created by the
 * python JSON serializer.
 */
function parse(result) {
  if (typeof result === 'string') {
    try {
      return JSON.parse(result)
    } catch (e) {
      try {
        return JSON.parse(result.replace(/\bNaN\b/g, '"NaN"'))
      } catch (e) {
        return result
      }
    }
  } else {
    return result
  }
}

export class ApiProviderComponent extends React.Component {
  static propTypes = {
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired,
    raiseError: PropTypes.func.isRequired,
    keycloak: PropTypes.object.isRequired,
    keycloakInitialized: PropTypes.bool
  }

  constructor(props) {
    super(props)
    this.onToken = this.onToken.bind(this)
  }

  onToken(token) {
    // console.log(token)
  }

  update() {
    const { keycloak } = this.props
    this.setState({api: this.createApi(keycloak)})
    if (keycloak.token) {
      keycloak.loadUserInfo()
        .success(user => {
          this.setState({user: user})
        })
        .error(error => {
          this.props.raiseError(error)
        })
    }
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (this.props.keycloakInitialized !== prevProps.keycloakInitialized) {
      this.update()
    }
  }

  createApi(keycloak) {
    const api = new Api(keycloak)
    return api
  }

  state = {
    apiv1: null
  }

  render() {
    const { children } = this.props
    return (
      <apiContext.Provider value={this.state}>
        {children}
      </apiContext.Provider>
    )
  }
}

class LoginRequiredUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    message: PropTypes.string
  }

  static styles = theme => ({
    root: {
      display: 'flex',
      alignItems: 'center',
      padding: theme.spacing(2),
      '& p': {
        marginRight: theme.spacing(2)
      }
    }
  })

  render() {
    const {classes, message} = this.props

    let loginMessage = ''
    if (message) {
      loginMessage = <Typography>
        {this.props.message}
      </Typography>
    }

    return (
      <div className={classes.root}>
        <div>
          {loginMessage}
        </div>
        <LoginLogout color="primary" />
      </div>
    )
  }
}

export function DisableOnLoading({children}) {
  const containerRef = useRef(null)
  const {apiv1} = useContext(apiContext)
  const handleLoading = useCallback((loading) => {
    const enable = loading ? 'none' : ''
    containerRef.current.style.pointerEvents = enable
    containerRef.current.style.userSelects = enable
  }, [])

  useEffect(() => {
    apiv1.onLoading(handleLoading)
    return () => {
      apiv1.removeOnLoading(handleLoading)
    }
  }, [apiv1, handleLoading])

  return <div ref={containerRef}>{children}</div>
}
DisableOnLoading.propTypes = {
  children: PropTypes.any.isRequired
}

export const ApiV1Provider = compose(withKeycloak, withErrors)(ApiProviderComponent)

const LoginRequired = withStyles(LoginRequiredUnstyled.styles)(LoginRequiredUnstyled)

const __reauthorize_trigger_changes = ['apiv1', 'calcId', 'uploadId', 'calc_id', 'upload_id']

class WithApiComponent extends React.Component {
  static propTypes = {
    raiseError: PropTypes.func.isRequired,
    loginRequired: PropTypes.bool,
    showErrorPage: PropTypes.bool,
    loginMessage: PropTypes.string,
    api: PropTypes.object,
    user: PropTypes.object,
    Component: PropTypes.any
  }

  state = {
    notAuthorized: false
  }

  constructor(props) {
    super(props)
    this.raiseError = this.raiseError.bind(this)
  }

  componentDidUpdate(prevProps) {
    if (__reauthorize_trigger_changes.find(key => this.props[key] !== prevProps[key])) {
      this.setState({notAuthorized: false})
    }
  }

  raiseError(error) {
    const { raiseError, showErrorPage } = this.props

    console.error(error)

    if (!showErrorPage) {
      raiseError(error)
    } else {
      if (error.name === 'NotAuthorized') {
        this.setState({notAuthorized: true})
      } else {
        raiseError(error)
      }
    }
  }

  render() {
    const { raiseError, loginRequired, loginMessage, Component, ...rest } = this.props
    const { apiv1, keycloak } = rest
    const { notAuthorized } = this.state
    if (notAuthorized) {
      if (keycloak.authenticated) {
        return (
          <div style={{marginTop: 24}}>
            <Typography variant="h6">Not Authorized</Typography>
            <Typography>
              You are not authorized to access this information. If someone send
              you a link to this data, ask the authors to make the data publicly available
              or share it with you.
            </Typography>
          </div>
        )
      } else {
        return (
          <LoginRequired message="You need to be logged in to access this information." />
        )
      }
    } else {
      if (apiv1) {
        if (keycloak.authenticated || !loginRequired) {
          return <Component {...rest} raiseError={this.raiseError} />
        } else {
          return <LoginRequired message={loginMessage} />
        }
      } else {
        return ''
      }
    }
  }
}

const WithKeycloakWithApiCompnent = withKeycloak(WithApiComponent)

export function withApi(loginRequired, showErrorPage, loginMessage) {
  return function(Component) {
    return withErrors(props => (
      <apiContext.Consumer>
        {apiContext => (
          <WithKeycloakWithApiCompnent
            loginRequired={loginRequired}
            loginMessage={loginMessage}
            showErrorPage={showErrorPage}
            Component={Component}
            {...props} {...apiContext}
          />
        )}
      </apiContext.Consumer>
    ))
  }
}
