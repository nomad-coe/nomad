/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import React, { useContext, useEffect, useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import { DoesNotExist } from '../api'
import { useDataStore } from '../DataStore'
import { useHistory } from 'react-router-dom'
import { getUrl } from '../nav/Routes'
import { apiBase } from '../../config'

export const uploadPageContext = React.createContext()

export function useUploadPageContext() {
  return useContext(uploadPageContext)
}

const UploadPageContext = React.memo(function UploadPageContext({uploadId, children}) {
  const dataStore = useDataStore()
  dataStore.resetIfNeeded(uploadId)
  const [uploadStoreObj, setUploadStoreObj] = useState(dataStore.getUpload(apiBase, uploadId))

  const onUploadStoreUpdated = useCallback((oldStoreObj, newStoreObj) => {
    setUploadStoreObj(newStoreObj)
  }, [setUploadStoreObj])

  useEffect(() => {
    return dataStore.subscribeToUpload(apiBase, uploadId, onUploadStoreUpdated, true, true)
  }, [dataStore, uploadId, onUploadStoreUpdated])

  const history = useHistory()

  useEffect(() => {
    if (uploadStoreObj.error && uploadStoreObj.error instanceof DoesNotExist && uploadStoreObj.deletionRequested) {
      history.push(getUrl('user/uploads', new URL(window.location.href).pathname))
    }
  })

  return <uploadPageContext.Provider value={uploadStoreObj}>
    {children}
  </uploadPageContext.Provider>
})
UploadPageContext.propTypes = {
  uploadId: PropTypes.string.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default UploadPageContext
