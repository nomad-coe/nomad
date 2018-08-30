import { UploadRequest } from '@navjobs/upload'
import { apiBase } from './config'

const networkError = () => {
  throw Error('Network related error, cannot reach API or object storage.')
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
              'Content-Type': 'application/gzip'
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
    if (this.proc.current_task_name !== this.proc.task_names[0]) {
      this.uploading = 100
    } else if (!created && this.uploading === null) {
      // if data came from server during a normal get (not create) and its still uploading
      // and the uploading is also not controlled locally then it ought to be a failure/abort
      this.proc.status = 'FAILURE'
      this.is_ready = true
      this.proc.errors = ['upload failed, probably aborted']
    }
  }

  update() {
    if (this.uploading !== null && this.uploading !== 100) {
      return new Promise(resolve => resolve(this))
    } else {
      return fetch(`${apiBase}/uploads/${this.upload_id}`)
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
      'Content-Type': 'application/json'
    }
  }
  return fetch(`${apiBase}/uploads`, fetchData)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
    .then(uploadJson => new Upload(uploadJson, true))
}

function getUploads() {
  return fetch(`${apiBase}/uploads`)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
    .then(uploadsJson => uploadsJson.map(uploadJson => new Upload(uploadJson)))
}

function archive(uploadHash, calcHash) {
  return fetch(`${apiBase}/archive/${uploadHash}/${calcHash}`)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function repo(uploadHash, calcHash) {
  return fetch(`${apiBase}/repo/${uploadHash}/${calcHash}`)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function repoAll(page, perPage) {
  return fetch(`${apiBase}/repo?page=${page}&per_page=${perPage}`)
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

function deleteUpload(uploadId) {
  return fetch(`${apiBase}/uploads/${uploadId}`, {method: 'DELETE'})
    .catch(networkError)
    .then(handleResponseErrors)
    .then(response => response.json())
}

const api = {
  createUpload: createUpload,
  deleteUpload: deleteUpload,
  getUploads: getUploads,
  archive: archive,
  repo: repo,
  repoAll: repoAll
}

export default api
