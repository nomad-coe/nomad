import { UploadRequest } from '@navjobs/upload'
import Swagger from 'swagger-client'
import { apiBase } from './config'

const auth_headers = {
  Authorization: 'Basic ' + btoa('sheldon.cooper@nomad-fairdi.tests.de:password')
}

const swaggerPromise = Swagger(`${apiBase}/swagger.json`, {
  authorizations: {
    // my_query_auth: new ApiKeyAuthorization('my-query', 'bar', 'query'),
    // my_header_auth: new ApiKeyAuthorization('My-Header', 'bar', 'header'),
    'HTTP Basic': {
      username: 'sheldon.cooper@nomad-fairdi.tests.de',
      password: 'password'
    }
    // cookie_: new CookieAuthorization('one=two')
  }
})

export class DoesNotExist extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'DoesNotExist'
  }
}

const handleApiError = (e) => {
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

const upload_to_gui_ids = {}
let gui_upload_id_counter = 0

class Upload {
  constructor(json) {
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
              ...auth_headers
            }
          },
          files: [file],
          progress: value => {
            this.uploading = value
          }
        }
      )
      if (uploadRequest.error) {
        handleApiError(uploadRequest.error)
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
        return swaggerPromise.then(client => client.apis.uploads.get_upload({
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

function createUpload(name) {
  return new Upload({
    name: name,
    tasks: ['uploading', 'extract', 'parse_all', 'cleanup'],
    current_task: 'uploading',
    uploading: 0,
    create_time: new Date()
  }, true)
}

async function getUploads() {
  const client = await swaggerPromise
  return client.apis.uploads.get_uploads()
    .catch(handleApiError)
    .then(response => response.body.map(uploadJson => {
      const upload = new Upload(uploadJson)
      upload.uploading = 100
      return upload
    }))
}

async function archive(uploadId, calcId) {
  const client = await swaggerPromise
  return client.apis.archive.get_archive_calc({
    upload_id: uploadId,
    calc_id: calcId
  })
    .catch(handleApiError)
    .then(response => response.body)
}

async function calcProcLog(uploadId, calcId) {
  const client = await swaggerPromise
  return client.apis.archive.get_archive_logs({
    upload_id: uploadId,
    calc_id: calcId
  })
    .catch(handleApiError)
    .then(response => response.text)
}

async function repo(uploadId, calcId) {
  const client = await swaggerPromise
  return client.apis.repo.get_repo_calc({
    upload_id: uploadId,
    calc_id: calcId
  })
    .catch(handleApiError)
    .then(response => response.body)
}

async function repoAll(page, perPage, owner) {
  const client = await swaggerPromise
  return client.apis.repo.get_calcs({
    page: page,
    per_page: perPage,
    ower: owner || 'all'
  })
    .catch(handleApiError)
    .then(response => response.body)
}

async function deleteUpload(uploadId) {
  const client = await swaggerPromise
  return client.apis.uploads.delete_upload({upload_id: uploadId})
    .catch(handleApiError)
    .then(response => response.body)
}

async function commitUpload(uploadId) {
  const client = await swaggerPromise
  return client.apis.uploads.exec_upload_command({
    upload_id: uploadId,
    payload: {
      command: 'commit'
    }
  })
    .catch(handleApiError)
    .then(response => response.body)
}

let cachedMetaInfo = null

async function getMetaInfo() {
  if (cachedMetaInfo) {
    return cachedMetaInfo
  } else {
    const loadMetaInfo = async(path) => {
      const client = await swaggerPromise
      return client.apis.archive.get_metainfo({metainfo_path: path})
        .catch(handleApiError)
        .then(response => response.body)
        .then(data => {
          if (!cachedMetaInfo) {
            cachedMetaInfo = {
              loadedDependencies: {}
            }
          }
          cachedMetaInfo.loadedDependencies[path] = true
          if (data.metaInfos) {
            data.metaInfos.forEach(info => {
              cachedMetaInfo[info.name] = info
              info.relativePath = path
            })
          }
          if (data.dependencies) {
            data.dependencies
              .filter(dep => cachedMetaInfo.loadedDependencies[dep.relativePath] !== true)
              .forEach(dep => {
                loadMetaInfo(dep.relativePath)
              })
          }
        })
    }
    await loadMetaInfo('all.nomadmetainfo.json')
    return cachedMetaInfo
  }
}

async function getUploadCommand() {
  const client = await swaggerPromise
  return client.apis.uploads.get_upload_command()
    .catch(handleApiError)
    .then(response => response.body.upload_command)
}

const api = {
  getUploadCommand: getUploadCommand,
  createUpload: createUpload,
  deleteUpload: deleteUpload,
  commitUpload: commitUpload,
  getUploads: getUploads,
  archive: archive,
  calcProcLog: calcProcLog,
  repo: repo,
  repoAll: repoAll,
  getMetaInfo: getMetaInfo
}

export default api
