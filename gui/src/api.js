import { apiBase } from './config'

const networkError = () => {
  throw Error('Network related error, cannot reach API or object storage.')
}

const handleResponseErrors = (response) => {
  if (!response.ok) {
    throw Error(`API/object storage error (${response.statusCode}): ${response.statusText}`)
  }
  return response
}

class Upload {
  constructor(json) {
    Object.assign(this, json)
  }

  uploadFile(file) {
    console.assert(this.presigned_url)
    return fetch(this.presigned_url, {
      method: 'PUT',
      headers: {
        'Content-Type': 'application/gzip'
      },
      body: file
    })
      .catch(networkError)
      .then(handleResponseErrors)
      .then(() => this)
  }

  update() {
    return fetch(`${apiBase}/uploads/${this.upload_id}`)
      .catch(networkError)
      .then(handleResponseErrors)
      .then(response => response.json())
      .then(uploadJson => {
        Object.assign(this, uploadJson)
        return this
      })
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
    .then(uploadJson => new Upload(uploadJson))
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
