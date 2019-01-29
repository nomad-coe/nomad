import React from 'react'
import PropTypes from 'prop-types'
import { withErrors } from './errors'
import { UploadRequest } from '@navjobs/upload'
import Swagger from 'swagger-client'
import { apiBase } from '../config'
import { Typography, withStyles } from '@material-ui/core'
import LoginLogout from './LoginLogout'

const ApiContext = React.createContext()

export class DoesNotExist extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'DoesNotExist'
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
    user = user || {}
    this.auth_headers = {
      'X-Token': user.token
    }
    this.swaggerPromise = Api.createSwaggerClient(user.token)

    this.handleApiError = this.handleApiError.bind(this)
  }

  handleApiError(e) {
    if (e.response) {
      const body = e.response.body
      const message = (body && body.message) ? body.message : e.response.statusText
      if (e.response.status === 404) {
        throw new DoesNotExist(message)
      } else {
        throw Error(`API error (${e.response.status}): ${message}`)
      }
    } else {
      throw Error('Network related error, cannot reach API: ' + e)
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
    const client = await this.swaggerPromise
    return client.apis.uploads.get_uploads()
      .catch(this.handleApiError)
      .then(response => response.body.map(uploadJson => {
        const upload = new Upload(uploadJson, this)
        upload.uploading = 100
        return upload
      }))
  }

  async archive(uploadId, calcId) {
    const client = await this.swaggerPromise
    return client.apis.archive.get_archive_calc({
      upload_id: uploadId,
      calc_id: calcId
    })
      .catch(this.handleApiError)
      .then(response => response.body)
  }

  async calcProcLog(uploadId, calcId) {
    const client = await this.swaggerPromise
    return client.apis.archive.get_archive_logs({
      upload_id: uploadId,
      calc_id: calcId
    })
      .catch(this.handleApiError)
      .then(response => response.text)
  }

  async repo(uploadId, calcId) {
    const client = await this.swaggerPromise
    return client.apis.repo.get_repo_calc({
      upload_id: uploadId,
      calc_id: calcId
    })
      .catch(this.handleApiError)
      .then(response => response.body)
  }

  async repoAll(page, perPage, owner) {
    const client = await this.swaggerPromise
    return client.apis.repo.get_calcs({
      page: page,
      per_page: perPage,
      ower: owner || 'all'
    })
      .catch(this.handleApiError)
      .then(response => response.body)
  }

  async deleteUpload(uploadId) {
    const client = await this.swaggerPromise
    return client.apis.uploads.delete_upload({upload_id: uploadId})
      .catch(this.handleApiError)
      .then(response => response.body)
  }

  async commitUpload(uploadId) {
    const client = await this.swaggerPromise
    return client.apis.uploads.exec_upload_command({
      upload_id: uploadId,
      payload: {
        command: 'commit'
      }
    })
      .catch(this.handleApiError)
      .then(response => response.body)
  }

  static async authenticate(userName, password) {
    const client = await this.swaggerPromise
    return client.apis.auth.get_token()
      .catch(error => {
        if (error.response.status !== 401) {
          this.handleApiError(error)
        }
      })
      .then(response => response !== undefined)
  }

  _cachedMetaInfo = null

  async getMetaInfo() {
    if (this._cachedMetaInfo) {
      return this._cachedMetaInfo
    } else {
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
      return this._cachedMetaInfo
    }
  }

  async getUploadCommand() {
    const client = await this.swaggerPromise
    return client.apis.uploads.get_upload_command()
      .catch(this.handleApiError)
      .then(response => response.body.upload_command)
  }
}

export class ApiProvider extends React.Component {
  static propTypes = {
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired
  }

  state = {
    api: new Api(),
    user: null,
    login: (userName, password, successCallback) => {
      Api.createSwaggerClient(userName, password)
        .catch(this.state.api.handleApiError)
        .then(client => {
          client.apis.auth.get_user()
            .catch(error => {
              if (error.response.status !== 401) {
                this.handleApiError(error)
              }
            })
            .then(response => {
              if (response) {
                const user = response.body
                this.setState({api: new Api(user), user: user})
                successCallback(true)
              } else {
                successCallback(false)
              }
            })
        })
    },
    logout: () => {
      this.setState({api: new Api(), user: null})
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
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      display: 'flex',
      alignItems: 'center',
      '& p': {
        marginRight: theme.spacing.unit * 2
      }
    }
  })

  render() {
    const {classes} = this.props
    return (
      <div className={classes.root}>
        <Typography>
          To upload data, you must have a nomad account and you must be logged in.
        </Typography>
        <LoginLogout variant="outlined" color="primary"/>
      </div>
    )
  }
}

const LoginRequired = withStyles(LoginRequiredUnstyled.styles)(LoginRequiredUnstyled)

export function withApi(loginRequired) {
  return function(Component) {
    function WithApiComponent(props) {
      return (
        <ApiContext.Consumer>
          {apiContext => (
            (apiContext.user || !loginRequired)
              ? <Component
                {...props} api={apiContext.api} user={apiContext.user}
                login={apiContext.login} logout={apiContext.logout} />
              : <LoginRequired />
          )}
        </ApiContext.Consumer>
      )
    }
    return withErrors(WithApiComponent)
  }
}
