import React, { useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { appBase } from '../../config'
import { App, H5GroveProvider } from '@h5web/app'
import { useApi } from '../api'
import { useErrors } from '../errors'

import '@h5web/lib/dist/styles.css'
import '@h5web/app/dist/styles.css'

const H5Web = (props) => {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [filepath, setFilepath] = useState(false)
  const {upload_id, filename} = props
  const authParams = {}
  if (api.keycloak) {
    authParams.token = api.keycloak.token
  }

  useEffect(() => {
    if (filename.includes('.yaml') || filename.includes('.yml')) {
      api.get(`/uploads/${upload_id}/raw/${filename}`)
        .then(data => {
          const fname = data.match(/output:[ \t]*(.*nxs)/i)?.[1]
          setFilepath(upload_id.substring(0, 2) + '/' + upload_id + '/raw/' + fname)
        })
        .catch(raiseError)
    } else {
      setFilepath(upload_id.substring(0, 2) + '/' + upload_id + '/raw/' + filename)
    }
  }, [filepath, api, filename, upload_id, raiseError])

  if (filepath) {
    return (
      <H5GroveProvider
        url={appBase + '/h5grove/'}
        filepath={filepath}
        axiosParams={{file: filepath, upload_id: upload_id, ...authParams}}
      >
        <App />
      </H5GroveProvider>
    )
  }
  return 'Loading...'
}
H5Web.propTypes = {
  upload_id: PropTypes.string.isRequired,
  filename: PropTypes.string.isRequired
}

export default H5Web
