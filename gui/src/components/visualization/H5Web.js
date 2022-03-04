import React from 'react'
import PropTypes from 'prop-types'
import { appBase } from '../../config'
import { App, H5GroveProvider } from '@h5web/app'
import { useApi } from '../api'

import '@h5web/app/dist/style-lib.css'
import '@h5web/app/dist/style.css'

const H5Web = ({upload_id, filename}) => {
  const {api} = useApi()
  const filepath = upload_id.substring(0, 2) + '/' + upload_id + '/raw/' + filename
  const authParams = {}
  if (api.keycloak) {
    authParams.token = api.keycloak.token
  }
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
H5Web.propTypes = {
  upload_id: PropTypes.string.isRequired,
  filename: PropTypes.string.isRequired
}

export default H5Web
