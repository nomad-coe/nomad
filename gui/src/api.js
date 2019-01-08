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
    this._assignFromJson(json)
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

  _assignFromJson(uploadJson) {
    Object.assign(this, uploadJson)
    if (this.calcs) {
      this.calcs.results.forEach(calc => {
        const archiveId = calc.archive_id.split('/')
        calc.upload_hash = archiveId[0]
        calc.calc_hash = archiveId[1]
      })
    }
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
            this._assignFromJson(uploadJson)
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

async function archive(uploadHash, calcHash) {
  const client = await swaggerPromise
  return client.apis.archive.get_archive_calc({
      upload_hash: uploadHash,
      calc_hash: calcHash
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.body)
}

async function calcProcLog(uploadHash, calcHash) {
  const client = await swaggerPromise
  console.log(uploadHash + calcHash)
  return client.apis.archive.get_archive_logs({
    upload_hash: uploadHash,
    calc_hash: calcHash
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

async function repo(uploadHash, calcHash) {
  const client = await swaggerPromise
  return client.apis.repo.get_repo_calc({
      upload_hash: uploadHash,
      calc_hash: calcHash
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

async function unstageUpload(uploadId) {
  const client = await swaggerPromise
  return client.apis.uploads.exec_upload_command({
      upload_id: uploadId,
      payload: {
        operation: 'unstage'
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
            cachedMetaInfo = {}
          }
          if (data.dependencies) {
            data.dependencies.forEach(dep => {
              loadMetaInfo(dep.relativePath)
            })
          }
          if (data.metaInfos) {
            data.metaInfos.forEach(info => {
              cachedMetaInfo[info.name] = info
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
  unstageUpload: unstageUpload,
  getUploads: getUploads,
  archive: archive,
  calcProcLog: calcProcLog,
  repo: repo,
  repoAll: repoAll,
  getMetaInfo: getMetaInfo
}

export default api
