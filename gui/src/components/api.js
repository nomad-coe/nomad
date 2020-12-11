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
import { UploadRequest } from '@navjobs/upload'
import Swagger from 'swagger-client'
import { apiBase } from '../config'
import { Typography, withStyles } from '@material-ui/core'
import LoginLogout from './LoginLogout'
import { compose } from 'recompose'
import { withKeycloak } from 'react-keycloak'
import * as searchQuantities from '../searchQuantities.json'

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

const upload_to_gui_ids = {}
let gui_upload_id_counter = 0

class Upload {
  constructor(json, api) {
    this.api = api

    // Cannot use upload_id as key in GUI, because uploads don't have an upload_id
    // before upload is completed
    if (json.upload_id) {
      // instance from the API
      this.gui_upload_id = upload_to_gui_ids[json.upload_id]
      if (this.gui_upload_id === undefined) {
        // never seen in the GUI, needs a GUI id
        this.gui_upload_id = gui_upload_id_counter++
        upload_to_gui_ids[json.upload_id] = this.gui_upload_id
      }
    } else {
      // new instance, not from the API
      this.gui_upload_id = gui_upload_id_counter++
    }
    Object.assign(this, json)
  }

  uploadFile(file) {
    const uploadFileWithProgress = async() => {
      const authHeaders = await this.api.authHeaders()
      let uploadRequest = await UploadRequest(
        {
          request: {
            url: `${apiBase}/uploads/?name=${this.name}`,
            method: 'PUT',
            headers: {
              'Content-Type': 'application/gzip',
              ...authHeaders
            }
          },
          files: [file],
          progress: value => {
            this.uploading = value
          }
        }
      )
      if (uploadRequest.error) {
        handleApiError(uploadRequest.response ? uploadRequest.response.message : uploadRequest.error)
      }
      if (uploadRequest.aborted) {
        throw Error('User abort')
      }
      this.uploading = 100
      this.upload_id = uploadRequest.response.upload_id
      upload_to_gui_ids[this.upload_id] = this.gui_upload_id
    }

    return uploadFileWithProgress()
      .then(() => this)
  }

  get(page, perPage, orderBy, order) {
    if (this.uploading !== null && this.uploading !== 100) {
      return new Promise(resolve => resolve(this))
    } else {
      if (this.upload_id) {
        return this.api.swagger().then(client => client.apis.uploads.get_upload({
          upload_id: this.upload_id,
          page: page || 1,
          per_page: perPage || 5,
          order_by: orderBy || 'mainfile',
          order: order || -1
        }))
          .catch(handleApiError)
          .then(response => response.body)
          .then(uploadJson => {
            Object.assign(this, uploadJson)
            return this
          })
      } else {
        return new Promise(resolve => resolve(this))
      }
    }
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
  swagger() {
    if (this.keycloak.token) {
      const self = this
      return new Promise((resolve, reject) => {
        self.keycloak.updateToken()
          .success(() => {
            self._swaggerClient
              .then(swaggerClient => {
                swaggerClient.authorizations = {
                  'OpenIDConnect Bearer Token': `Bearer ${self.keycloak.token}`
                }
                resolve(swaggerClient)
              })
              .catch(() => {
                reject(new ApiError())
              })
          })
          .error(() => {
            reject(new ApiError())
          })
      })
    } else {
      const self = this
      return new Promise((resolve, reject) => {
        self._swaggerClient
          .then(swaggerClient => {
            swaggerClient.authorizations = {}
            resolve(swaggerClient)
          })
          .catch(() => {
            reject(new ApiError())
          })
      })
    }
  }

  authHeaders() {
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

  constructor(keycloak) {
    this._swaggerClient = Swagger(`${apiBase}/swagger.json`)
    this.keycloak = keycloak

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

  createUpload(name) {
    const upload = new Upload({
      upload_id: Api.uploadIds++,
      name: name,
      tasks: ['uploading', 'extract', 'parse_all', 'cleanup'],
      current_task: 'uploading',
      uploading: 0,
      create_time: new Date()
    }, this)

    return upload
  }

  onLoading(handler) {
    this.loadingHandler = [...this.loadingHandler, handler]
  }

  removeOnLoading(handler) {
    this.loadingHandler = this.loadingHandler.filter(item => item !== handler)
  }

  async getUploads(state, page, perPage) {
    state = state || 'all'
    page = page || 1
    perPage = perPage || 10

    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.get_uploads({state: state, page: page, per_page: perPage}))
      .catch(handleApiError)
      .then(response => ({
        ...response.body,
        results: response.body.results.map(uploadJson => {
          const upload = new Upload(uploadJson, this)
          upload.uploading = 100
          return upload
        })
      }))
      .finally(this.onFinishLoading)
  }

  async getUnpublishedUploads() {
    return this.getUploads('unpublished', 1, 1000)
  }

  async getPublishedUploads(page, perPage) {
    return this.getUploads('published', 1, 10)
  }

  async archive(uploadId, calcId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.archive.get_archive_calc({
        upload_id: uploadId,
        calc_id: calcId
      }))
      .catch(handleApiError)
      .then(response => {
        const result = response.body || response.text || response.data
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
      })
      .finally(this.onFinishLoading)
  }

  async encyclopediaBasic(materialId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.encyclopedia.get_material({
        material_id: materialId
      }))
      .catch(handleApiError)
      .then(response => {
        const result = response.body || response.text || response.data
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
      })
      .finally(this.onFinishLoading)
  }

  async encyclopediaCalculations(materialId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.encyclopedia.get_calculations({
        material_id: materialId
      }))
      .catch(handleApiError)
      .then(response => {
        const result = response.body || response.text || response.data
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
      })
      .finally(this.onFinishLoading)
  }

  async encyclopediaCalculation(materialId, calcId, payload) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.encyclopedia.get_calculation({
        material_id: materialId,
        calc_id: calcId,
        payload: payload
      }))
      .catch(handleApiError)
      .then(response => {
        const result = response.body || response.text || response.data
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
      })
      .finally(this.onFinishLoading)
  }

  async calcProcLog(uploadId, calcId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.archive.get_archive_logs({
        upload_id: uploadId,
        calc_id: calcId
      }))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async getRawFileListFromCalc(uploadId, calcId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.raw.get_file_list_from_calc({
        upload_id: uploadId,
        calc_id: calcId,
        path: null
      }))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async getRawFile(uploadId, calcId, path, kwargs) {
    this.onStartLoading()
    const length = (kwargs && kwargs.length) || 4096
    return this.swagger()
      .then(client => client.apis.raw.get_file_from_calc({
        upload_id: uploadId,
        calc_id: calcId,
        path: path,
        decompress: true,
        ...(kwargs || {}),
        length: length
      }))
      .catch(handleApiError)
      .then(response => {
        /* global Blob */
        /* eslint no-undef: "error" */
        if (response.data instanceof Blob) {
          if (response.data.type.endsWith('empty')) {
            return {
              contents: '',
              hasMore: false,
              mimeType: 'plain/text'
            }
          }
          return {
            contents: null,
            hasMore: false,
            mimeType: response.data.type
          }
        }
        return {
          contents: response.data,
          hasMore: response.data.length === length,
          mimeType: 'plain/text'
        }
      })
      .finally(this.onFinishLoading)
  }

  async repo(uploadId, calcId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.repo.get_repo_calc({
        upload_id: uploadId,
        calc_id: calcId
      }))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async edit(edit) {
    // We do not call the start and finish loading callbacks, because this one is
    // only used in the background.

    // repair the query, the API will only access correct use of lists for many
    // quantities
    Object.keys(edit.query).forEach(quantity => {
      if (searchQuantities[quantity] && searchQuantities[quantity].many) {
        if (!Array.isArray(edit.query[quantity])) {
          edit.query[quantity] = edit.query[quantity].split(',')
        }
      }
    })
    return this.swagger()
      .then(client => client.apis.repo.edit_repo({payload: edit}))
      .catch(handleApiError)
      .then(response => response.body)
  }

  async resolvePid(pid) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.repo.resolve_pid({pid: pid}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async resolveDoi(doi) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.datasets.resolve_doi({doi: doi}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async search(search) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.repo.search({
        exclude: ['atoms', 'only_atoms', 'dft.files', 'dft.quantities', 'dft.optimade', 'dft.labels', 'dft.geometries'],
        ...search}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async getDatasets(prefix) {
    // no loading indicator, because this is only used in the background of the edit dialog
    return this.swagger()
      .then(client => client.apis.datasets.list_datasets({prefix: prefix}))
      .catch(handleApiError)
      .then(response => response.body)
  }

  async assignDatasetDOI(datasetName) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.datasets.assign_doi({name: datasetName}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async deleteDataset(datasetName) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.datasets.delete_dataset({name: datasetName}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async getUsers(query) {
    // no loading indicator, because this is only used in the background of the edit dialog
    return this.swagger()
      .then(client => client.apis.auth.get_users({query: query}))
      .catch(handleApiError)
      .then(response => response.body)
  }

  async inviteUser(user) {
    return this.swagger()
      .then(client => client.apis.auth.invite_user({payload: user}))
      .catch(handleApiError)
      .then(response => response.body)
  }

  async quantities_search(search) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.repo.quantities_search(search))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async quantity_search(search) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.repo.quantity_search(search))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async suggestions_search(quantity, search, include, size, noLoadingIndicator) {
    if (!noLoadingIndicator) {
      this.onStartLoading()
    }
    return this.swagger()
      .then(client => client.apis.repo.suggestions_search({
        size: size || 20,
        include: include,
        quantity: quantity,
        ...search
      }))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(() => {
        if (!noLoadingIndicator) {
          this.onFinishLoading()
        }
      })
  }

  async deleteUpload(uploadId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.delete_upload({upload_id: uploadId}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async publishUpload(uploadId, embargoLength) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.exec_upload_operation({
        upload_id: uploadId,
        payload: {
          operation: 'publish',
          metadata: {
            with_embargo: embargoLength > 0,
            embargo_length: embargoLength
          }
        }
      }))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async publishUploadToCentralNomad(uploadId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.exec_upload_operation({
        upload_id: uploadId,
        payload: {
          operation: 'publish-to-central-nomad'
        }
      }))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async getSignatureToken() {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.auth.get_auth())
      .catch(handleApiError)
      .then(response => response.body.signature_token)
      .finally(this.onFinishLoading)
  }

  _cachedInfo = null

  async getInfo() {
    if (!this._cachedInfo) {
      this.onStartLoading()
      this._cachedInfo = await this.swagger()
        .then(client => {
          return client.apis.info.get_info()
            .then(response => response.body)
            .catch(handleApiError)
        })
        .finally(this.onFinishLoading)
    }
    return this._cachedInfo
  }

  async getUploadCommand() {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.get_upload_command())
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
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
    api.getInfo()
      .catch(handleApiError)
      .then(info => {
        if (info.parsers) {
          info.parsers.sort()
        }
        this.setState({info: info})
      })
      .catch(error => {
        this.props.raiseError(error)
      })

    return api
  }

  state = {
    api: null,
    info: null
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
  const {api} = useContext(apiContext)
  const handleLoading = useCallback((loading) => {
    const enable = loading ? 'none' : ''
    containerRef.current.style.pointerEvents = enable
    containerRef.current.style.userSelects = enable
  }, [])

  useEffect(() => {
    api.onLoading(handleLoading)
    return () => {
      api.removeOnLoading(handleLoading)
    }
  }, [api, handleLoading])

  return <div ref={containerRef}>{children}</div>
}
DisableOnLoading.propTypes = {
  children: PropTypes.any.isRequired
}

export const ApiProvider = compose(withKeycloak, withErrors)(ApiProviderComponent)

const LoginRequired = withStyles(LoginRequiredUnstyled.styles)(LoginRequiredUnstyled)

const __reauthorize_trigger_changes = ['api', 'calcId', 'uploadId', 'calc_id', 'upload_id']

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
    const { api, keycloak } = rest
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
      if (api) {
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
