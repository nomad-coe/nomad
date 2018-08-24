import { apiBase } from './config'

class Upload {
  constructor(json) {
    console.debug('Created local upload for ' + json.upload_id)
    Object.assign(this, json)
  }

  uploadFile(file) {
    console.assert(this.presigned_url)
    console.debug(`Upload ${file} to ${this.presigned_url}.`)
    return fetch(this.presigned_url, {
      method: 'PUT',
      headers: {
        'Content-Type': 'application/gzip'
      },
      body: file
    }).then(() => {
      console.log(`Uploaded ${file} to ${this.upload_id}.`)
      return this
    }).catch(error => {
      console.error(`Could not upload ${file} to ${this.presigned_url}: ${error}.`)
      return this
    })
  }

  update() {
    return fetch(`${apiBase}/uploads/${this.upload_id}`)
      .then(response => response.json())
      .then(uploadJson => {
        Object.assign(this, uploadJson)
        return this
      })
  }
}

function createUpload(name) {
  console.debug('Request new upload.')
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
    .then(response => response.json())
    .then(uploadJson => new Upload(uploadJson))
}

function getUploads() {
  return fetch(`${apiBase}/uploads`)
    .then(response => response.json())
    .then(uploadsJson => uploadsJson.map(uploadJson => new Upload(uploadJson)))
}

function archive(uploadHash, calcHash) {
  return fetch(`${apiBase}/archive/${uploadHash}/${calcHash}`)
    .then(response => response.json())
}

const api = {
  createUpload: createUpload,
  getUploads: getUploads,
  archive: archive
};

export default api;