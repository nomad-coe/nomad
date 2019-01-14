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

const networkError = (e) => {
  console.log(e)
  throw Error('Network related error, cannot reach API: ' + e)
}

const handleJsonErrors = (e) => {
  console.log(e)
  throw Error('API return unexpected data format.')
}

const handleResponseErrors = (response) => {
  if (!response.ok) {
    return response.json()
      .catch(() => {
        throw Error(`API error (${response.status}): ${response.statusText}`)
      }).then(data => {
        throw Error(`API error (${response.status}): ${data.message}`)
      })
  }
  return response
}

class Upload {
  constructor(json) {
    this.uploading = 0
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
        networkError(uploadRequest.error)
      }
      if (uploadRequest.aborted) {
        throw Error('User abort')
      }
      this.uploading = 100
      this.upload_id = uploadRequest.response.upload_id
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
            page: page || 1,
            per_page: perPage || 5,
            order_by: orderBy || 'mainfile',
            order: order || -1
          }))
          .catch(networkError)
          .then(handleResponseErrors)
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
    tasks: ['uploading'],
    current_task: 'uploading',
    uploading: 0,
    create_time: new Date()
  }, true)
}

async function getUploads() {
  const client = await swaggerPromise
  return client.apis.uploads.get_uploads()
    .catch(networkError)
    .then(handleResponseErrors)
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
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.body)
}

async function calcProcLog(uploadId, calcId) {
  const client = await swaggerPromise
  return client.apis.archive.get_archive_logs({
    upload_id: uploadId,
    calc_id: calcId
  })
    .catch(networkError)
    .then(response => {
      if (!response.ok) {
        if (response.status === 404) {
          return ''
        } else {
          return handleResponseErrors(response)
        }
      } else {
        return response.text
      }
    })
}

async function repo(uploadId, calcId) {
  const client = await swaggerPromise
  return client.apis.repo.get_repo_calc({
      upload_id: uploadId,
      calc_id: calcId
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.body)
}

async function repoAll(page, perPage, owner) {
  const client = await swaggerPromise
  return client.apis.repo.get_calcs({
      page: page,
      per_page: perPage,
      ower: owner || 'all'
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.body)
}

async function deleteUpload(uploadId) {
  const client = await swaggerPromise
  return client.apis.uploads.delete_upload({upload_id: uploadId})
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.body)
}

async function commitUpload(uploadId) {
  const client = await swaggerPromise
  return client.apis.uploads.exec_upload_command({
      upload_id: uploadId,
      payload: {
        operation: 'commit'
      }
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.body)
}

let cachedMetaInfo = null

async function getMetaInfo() {
  if (cachedMetaInfo) {
    return cachedMetaInfo
  } else {
    const loadMetaInfo = async(path) => {
      return fetch(`${apiBase}/archive/metainfo/${path}`)
        .catch(networkError)
        .then(handleResponseErrors)
        .then(response => response.json())
        .catch(handleJsonErrors)
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
    .catch(networkError)
    .then(handleResponseErrors)
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
