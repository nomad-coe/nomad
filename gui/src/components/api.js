import React from 'react'
import PropTypes from 'prop-types'
import { withErrors } from './errors'
import { UploadRequest } from '@navjobs/upload'
import Swagger from 'swagger-client'
import { apiBase } from '../config'
import { Typography, withStyles, Link } from '@material-ui/core'
import LoginLogout from './LoginLogout'
import { compose } from 'recompose'
import MetaInfoRepository from './MetaInfoRepository'
import { withKeycloak } from 'react-keycloak'

const ApiContext = React.createContext()

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
    const message = (body && body.message) ? body.message : e.response.statusText
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
    this.onStartLoading = () => null
    this.onFinishLoading = () => null

    this._swaggerClient = Swagger(`${apiBase}/swagger.json`)
    this.keycloak = keycloak

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

  async calcProcLog(uploadId, calcId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.archive.get_archive_logs({
        upload_id: uploadId,
        calc_id: calcId
      }))
      .catch(handleApiError)
      .then(response => response.text)
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

  async getRawFile(uploadId, path, kwargs) {
    this.onStartLoading()
    const length = (kwargs && kwargs.length) || 4096
    return this.swagger()
      .then(client => client.apis.raw.get({
        upload_id: uploadId,
        path: path,
        decompress: true,
        ...(kwargs || {}),
        length: length
      }))
      .catch(handleApiError)
      .then(response => {
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
    // this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.repo.edit_repo({payload: edit}))
      .catch(handleApiError)
      .then(response => response.body)
      // .finally(this.onFinishLoading)
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
      .then(client => client.apis.repo.search(search))
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

  async deleteUpload(uploadId) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.delete_upload({upload_id: uploadId}))
      .catch(handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async publishUpload(uploadId, withEmbargo) {
    this.onStartLoading()
    return this.swagger()
      .then(client => client.apis.uploads.exec_upload_operation({
        upload_id: uploadId,
        payload: {
          operation: 'publish',
          metadata: {
            with_embargo: withEmbargo
          }
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

  _metaInfoRepositories = {}

  async getMetaInfo(pkg) {
    pkg = pkg || 'all.nomadmetainfo.json'

    const metaInfoRepository = this._metaInfoRepositories[pkg]

    if (metaInfoRepository) {
      return metaInfoRepository
    } else {
      this.onStartLoading()
      try {
        const loadMetaInfo = async(path) => {
          return this.swagger()
            .then(client => client.apis.archive.get_metainfo({metainfo_package_name: path}))
            .catch(handleApiError)
            .then(response => response.body)
        }
        const metaInfo = await loadMetaInfo(pkg)
        const metaInfoRepository = new MetaInfoRepository(metaInfo)
        this._metaInfoRepositories[pkg] = metaInfoRepository

        return metaInfoRepository
      } finally {
        this.onFinishLoading()
      }
    }
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
    api.onStartLoading = (name) => {
      this.setState(state => ({loading: state.loading + 1}))
    }
    api.onFinishLoading = (name) => {
      this.setState(state => ({loading: Math.max(0, state.loading - 1)}))
    }

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
    info: null,
    loading: 0
  }

  render() {
    const { children } = this.props
    return (
      <ApiContext.Provider value={this.state}>
        {children}
      </ApiContext.Provider>
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
      padding: theme.spacing.unit * 2,
      '& p': {
        marginRight: theme.spacing.unit * 2
      }
    }
  })

  render() {
    const {classes, message} = this.props

    let loginMessage = ''
    if (message) {
      loginMessage = <Typography>
        {this.props.message} If you do not have a NOMAD Repository account, register <Link href='http://nomad-coe.eu:8080/NomadRepository-1.1/register/'>here</Link>.
      </Typography>
    }

    return (
      <div className={classes.root}>
        {loginMessage}
        <LoginLogout variant="outlined" color="primary" />
      </div>
    )
  }
}

export function DisableOnLoading(props) {
  return (
    <ApiContext.Consumer>
      {apiContext => (
        <div style={apiContext.loading ? { pointerEvents: 'none', userSelects: 'none' } : {}}>
          {props.children}
        </div>
      )}
    </ApiContext.Consumer>
  )
}
DisableOnLoading.propTypes = {
  children: PropTypes.any.isRequired
}

export const ApiProvider = compose(withKeycloak, withErrors)(ApiProviderComponent)

const LoginRequired = withStyles(LoginRequiredUnstyled.styles)(LoginRequiredUnstyled)

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
    if (prevProps.api !== this.props.api) {
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
          <div>
            <Typography variant="h6">Not Authorized</Typography>
            <Typography>
              You are not authorized to access this information. If someone send
              you this link, ask him to make his data publicly available or share
              it with you.
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
      <ApiContext.Consumer>
        {apiContext => (
          <WithKeycloakWithApiCompnent
            loginRequired={loginRequired}
            loginMessage={loginMessage}
            showErrorPage={showErrorPage}
            Component={Component}
            {...props} {...apiContext}
          />
        )}
      </ApiContext.Consumer>
    ))
  }
}
