import { UploadRequest } from '@navjobs/upload'
import { apiBase, appStaticBase } from './config'

const auth_headers = {
  Authorization: 'Basic ' + btoa('me@gmail.com:nomad')
}

const networkError = () => {
  throw Error('Network related error, cannot reach API or object storage.')
}

const handleJsonErrors = () => {
  throw Error('Server return unexpected data format.')
}

const handleResponseErrors = (response) => {
  if (!response.ok) {
    return response.json()
      .catch(() => {
        throw Error(`API/object storage error (${response.status}): ${response.statusText}`)
      }).then(data => {
        throw Error(`API/object storage error (${response.status}): ${data.message}`)
      })
  }
  return response
}

class Upload {
  constructor(json, created) {
    this.uploading = null
    this._assignFromJson(json, created)
  }

  uploadFile(file) {
    console.assert(this.presigned_url)
    this.uploading = 0

    const uploadFileWithProgress = async() => {
      let { error, aborted } = await UploadRequest(
        {
          request: {
            url: this.presigned_url,
            method: 'PUT',
            headers: {
              'Content-Type': 'application/gzip',
              ...auth_headers
            }
          },
          files: [file],
          progress: value => {
            console.log(value)
            this.uploading = value
          }
        }
      )
      if (error) {
        networkError(error)
      }
      if (aborted) {
        throw Error('User abort')
      }
    }

    return uploadFileWithProgress()
      .then(() => this)
  }

  _assignFromJson(uploadJson, created) {
    Object.assign(this, uploadJson)
    if (this.current_task !== this.tasks[0]) {
      this.uploading = 100
      this.waiting = false
    } else if (!created && this.uploading === null) {
      // if data came from server during a normal get (not create) and its still uploading
      // and the uploading is also not controlled locally then it ought to be a manual upload
      this.waiting = true
    }
  }

  get(page, perPage, orderBy, order) {
    if (!page) page = 1
    if (!perPage) perPage = 5
    if (!orderBy) orderBy = 'mainfile'
    if (!order) order = 'desc'

    order = order === 'desc' ? -1 : 1

    if (this.uploading !== null && this.uploading !== 100) {
      return new Promise(resolve => resolve(this))
    } else {
      const qparams = `page=${page}&per_page=${perPage}&order_by=${orderBy}&order=${order}`
      return fetch(
        `${apiBase}/uploads/${this.upload_id}?${qparams}`,
        {
          method: 'GET',
          headers: auth_headers
        })
        .catch(networkError)
        .then(handleResponseErrors)
        .then(response => response.json())
        .then(uploadJson => {
          this._assignFromJson(uploadJson)
          return this
        })
    }
  }
}

function createUpload(name) {
  const fetchData = {
    method: 'POST',
    body: JSON.stringify({
      name: name
    }),
    headers: {
      'Content-Type': 'application/json',
      ...auth_headers
    }
  }
  return fetch(`${apiBase}/uploads`, fetchData)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
    .then(uploadJson => new Upload(uploadJson, true))
}

function getUploads() {
  return fetch(
    `${apiBase}/uploads`,
    {
      method: 'GET',
      headers: auth_headers
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
    .then(uploadsJson => uploadsJson.map(uploadJson => new Upload(uploadJson)))
}

function archive(uploadHash, calcHash) {
  return fetch(archiveUrl(uploadHash, calcHash))
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function archiveUrl(uploadHash, calcHash) {
  return `${apiBase}/archive/${uploadHash}/${calcHash}`
}

function repo(uploadHash, calcHash) {
  return fetch(`${apiBase}/repo/${uploadHash}/${calcHash}`)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function repoAll(page, perPage, owner) {
  return fetch(
    `${apiBase}/repo?page=${page}&per_page=${perPage}&owner=${owner || 'all'}`,
    {
      method: 'GET',
      headers: auth_headers
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function deleteUpload(uploadId) {
  return fetch(
    `${apiBase}/uploads/${uploadId}`,
    {
      method: 'DELETE',
      headers: auth_headers
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function unstageUpload(uploadId) {
  return fetch(
    `${apiBase}/uploads/${uploadId}`,
    {
      method: 'POST',
      body: JSON.stringify({
        operation: 'unstage'
      }),
      headers: {
        'Content-Type': 'application/json',
        ...auth_headers
      }
    })
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

let cachedMetaInfo = null

async function getMetaInfo() {
  if (cachedMetaInfo) {
    return cachedMetaInfo
  } else {
    const loadMetaInfo = async(path) => {
      return fetch(`${appStaticBase}/metainfo/meta_info/nomad_meta_info/${path}`)
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

const api = {
  createUpload: createUpload,
  deleteUpload: deleteUpload,
  unstageUpload: unstageUpload,
  getUploads: getUploads,
  archive: archive,
  archiveUrl: archiveUrl,
  repo: repo,
  repoAll: repoAll,
  getMetaInfo: getMetaInfo
}

export default api
