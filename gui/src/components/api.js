import React from 'react'
import PropTypes, { instanceOf } from 'prop-types'
import { withErrors } from './errors'
import { UploadRequest } from '@navjobs/upload'
import Swagger from 'swagger-client'
import { apiBase } from '../config'
import { Typography, withStyles } from '@material-ui/core'
import LoginLogout from './LoginLogout'
import { Cookies, withCookies } from 'react-cookie'
import { compose } from 'recompose'

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

const upload_to_gui_ids = {}
let gui_upload_id_counter = 0

class Upload {
  constructor(json, api) {
    this.api = api
    this.handleApiError = api.handleApiError.bind(api)

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
      let uploadRequest = await UploadRequest(
        {
          request: {
            url: `${apiBase}/uploads/?name=${this.name}`,
            method: 'PUT',
            headers: {
              'Content-Type': 'application/gzip',
              ...this.api.auth_headers
            }
          },
          files: [file],
          progress: value => {
            this.uploading = value
          }
        }
      )
      if (uploadRequest.error) {
        this.handleApiError(uploadRequest.error)
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
        return this.api.swaggerPromise.then(client => client.apis.uploads.get_upload({
          upload_id: this.upload_id,
          page: page || 1,
          per_page: perPage || 5,
          order_by: orderBy || 'mainfile',
          order: order || -1
        }))
          .catch(this.handleApiError)
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

class Api {
  static async createSwaggerClient(userNameToken, password) {
    let data
    if (userNameToken) {
      let auth = {
        'X-Token': userNameToken
      }
      if (password) {
        auth = {
          'HTTP Basic': {
            username: userNameToken,
            password: password
          }
        }
      }
      data = {authorizations: auth}
    }

    return Swagger(`${apiBase}/swagger.json`, data)
  }

  constructor(user) {
    this.onStartLoading = () => null
    this.onFinishLoading = () => null

    this.handleApiError = this.handleApiError.bind(this)

    user = user || {}
    this.auth_headers = {
      'X-Token': user.token
    }
    this.swaggerPromise = Api.createSwaggerClient(user.token)
  }

  handleApiError(e) {
    if (e.response) {
      const body = e.response.body
      const message = (body && body.message) ? body.message : e.response.statusText
      if (e.response.status === 404) {
        throw new DoesNotExist(message)
      } else if (e.response.status === 401) {
        throw new NotAuthorized(message)
      } else {
        throw Error(`API error (${e.response.status}): ${message}`)
      }
    } else {
      throw Error('Network related error, cannot reach API')
    }
  }

  createUpload(name) {
    return new Upload({
      name: name,
      tasks: ['uploading', 'extract', 'parse_all', 'cleanup'],
      current_task: 'uploading',
      uploading: 0,
      create_time: new Date()
    }, this)
  }

  async getUploads() {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.uploads.get_uploads()
      .catch(this.handleApiError)
      .then(response => response.body.map(uploadJson => {
        const upload = new Upload(uploadJson, this)
        upload.uploading = 100
        return upload
      }))
      .finally(this.onFinishLoading)
  }

  async archive(uploadId, calcId) {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.archive.get_archive_calc({
      upload_id: uploadId,
      calc_id: calcId
    })
      .catch(this.handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async calcProcLog(uploadId, calcId) {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.archive.get_archive_logs({
      upload_id: uploadId,
      calc_id: calcId
    })
      .catch(this.handleApiError)
      .then(response => response.text)
      .finally(this.onFinishLoading)
  }

  async repo(uploadId, calcId) {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.repo.get_repo_calc({
      upload_id: uploadId,
      calc_id: calcId
    })
      .catch(this.handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async search(search) {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.repo.search(search)
      .catch(this.handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async deleteUpload(uploadId) {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.uploads.delete_upload({upload_id: uploadId})
      .catch(this.handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async publishUpload(uploadId, withEmbargo) {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.uploads.exec_upload_operation({
      upload_id: uploadId,
      payload: {
        operation: 'publish',
        metadata: {
          with_embargo: withEmbargo
        }
      }
    })
      .catch(this.handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  async getSignatureToken() {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.auth.get_token()
      .catch(this.handleApiError)
      .then(response => response.body)
      .finally(this.onFinishLoading)
  }

  _cachedMetaInfo = null

  async getMetaInfo() {
    if (this._cachedMetaInfo) {
      return this._cachedMetaInfo
    } else {
      this.onStartLoading()
      const loadMetaInfo = async(path) => {
        const client = await this.swaggerPromise
        return client.apis.archive.get_metainfo({metainfo_path: path})
          .catch(this.handleApiError)
          .then(response => response.body)
          .then(data => {
            if (!this._cachedMetaInfo) {
              this._cachedMetaInfo = {
                loadedDependencies: {}
              }
            }
            this._cachedMetaInfo.loadedDependencies[path] = true
            if (data.metaInfos) {
              data.metaInfos.forEach(info => {
                this._cachedMetaInfo[info.name] = info
                info.relativePath = path
              })
            }
            if (data.dependencies) {
              data.dependencies
                .filter(dep => this._cachedMetaInfo.loadedDependencies[dep.relativePath] !== true)
                .forEach(dep => {
                  loadMetaInfo(dep.relativePath)
                })
            }
          })
      }
      await loadMetaInfo('all.nomadmetainfo.json')
      this.onFinishLoading()
      return this._cachedMetaInfo
    }
  }

  _cachedInfo = null

  async getInfo() {
    if (!this._cachedInfo) {
      this.onStartLoading()
      this._cachedInfo = this.swaggerPromise
        .then(client => {
          return client.apis.info.get_info()
            .catch(this.handleApiError)
            .then(response => {
              this.onFinishLoading()
              return response.body
            })
        })
    }
    return this._cachedInfo
  }

  async getUploadCommand() {
    this.onStartLoading()
    const client = await this.swaggerPromise
    return client.apis.uploads.get_upload_command()
      .catch(this.handleApiError)
      .then(response => response.body.upload_command)
      .finally(this.onFinishLoading)
  }
}

export class ApiProviderComponent extends React.Component {
  static propTypes = {
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired,
    cookies: instanceOf(Cookies).isRequired,
    raiseError: PropTypes.func.isRequired
  }

  componentDidMount(props) {
    const token = this.props.cookies.get('token')
    if (token && token !== 'undefined') {
      this.state.login(token)
    } else {
      this.setState({api: this.createApi()})
    }
  }

  createApi(user) {
    const api = new Api(user)
    api.onStartLoading = () => this.setState({loading: this.state.loading + 1})
    api.onFinishLoading = () => this.setState({loading: Math.max(0, this.state.loading - 1)})
    return api
  }

  state = {
    api: null,
    user: null,
    isLoggingIn: false,
    loading: 0,
    login: (userNameToken, password, successCallback) => {
      this.setState({isLoggingIn: true})
      successCallback = successCallback || (() => true)
      Api.createSwaggerClient(userNameToken, password)
        .catch((error) => {
          this.setState({api: this.createApi(), isLoggingIn: false, user: null})
          this.props.raiseError(error)
        })
        .then(client => {
          client.apis.auth.get_user()
            .catch(error => {
              if (error.response.status !== 401) {
                try {
                  this.props.raiseError(error)
                } catch (e) {
                  this.setState({api: this.createApi(), isLoggingIn: false, user: null})
                  this.props.raiseError(error)
                }
              }
            })
            .then(response => {
              if (response) {
                const user = response.body
                this.setState({api: this.createApi(user), isLoggingIn: false, user: user})
                this.props.cookies.set('token', user.token)
                successCallback(true)
              } else {
                this.setState({api: this.createApi(), isLoggingIn: false, user: null})
                successCallback(false)
              }
            })
        })
        .catch(error => {
          this.setState({api: this.createApi(), isLoggingIn: false, user: null})
          this.props.raiseError(error)
        })
    },
    logout: () => {
      this.setState({api: this.createApi(), user: null})
      this.props.cookies.set('token', undefined)
    }
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
    message: PropTypes.string,
    isLoggingIn: PropTypes.bool,
    onLoggedIn: PropTypes.func
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
    const {classes, isLoggingIn, onLoggedIn} = this.props
    return (
      <div className={classes.root}>
        <Typography>
          {this.props.message || ''}
        </Typography>
        <LoginLogout onLoggedIn={onLoggedIn} variant="outlined" color="primary" isLoggingIn={isLoggingIn}/>
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

export const ApiProvider = compose(withCookies, withErrors)(ApiProviderComponent)

const LoginRequired = withStyles(LoginRequiredUnstyled.styles)(LoginRequiredUnstyled)

class WithApiComponent extends React.Component {
  static propTypes = {
    raiseError: PropTypes.func.isRequired,
    loginRequired: PropTypes.bool,
    showErrorPage: PropTypes.bool,
    loginMessage: PropTypes.string,
    api: PropTypes.object,
    user: PropTypes.object,
    isLoggingIn: PropTypes.bool,
    Component: PropTypes.any
  }

  state = {
    notAuthorized: false,
    notFound: false
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

    if (!showErrorPage) {
      raiseError(error)
    } else {
      if (error.name === 'NotAuthorized') {
        this.setState({notAuthorized: true})
      } else if (error.name === 'DoesNotExist') {
        this.setState({notFound: true})
      } else {
        raiseError(error)
      }
    }
  }

  render() {
    const { raiseError, loginRequired, loginMessage, Component, ...rest } = this.props
    const { api, user, isLoggingIn } = rest
    const { notAuthorized, notFound } = this.state
    if (notAuthorized) {
      if (user) {
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
          <LoginRequired
            message="You need to be logged in to access this information."
            isLoggingIn={isLoggingIn}
            onLoggedIn={() => this.setState({notAuthorized: false})}
          />
        )
      }
    } else if (notFound) {
      return <div>
        <Typography variant="h6">Not Found</Typography>
        <Typography>
        The information that you are trying to access does not exists.
        </Typography>
      </div>
    } else {
      if (api) {
        if (user || !loginRequired) {
          return <Component {...rest} raiseError={this.raiseError} />
        } else {
          return <LoginRequired message={loginMessage} isLoggingIn={isLoggingIn} />
        }
      } else {
        return ''
      }
    }
  }
}

export function withApi(loginRequired, showErrorPage, loginMessage) {
  return function(Component) {
    return withErrors(props => (
      <ApiContext.Consumer>
        {apiContext => (
          <WithApiComponent
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
