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

import React, { useCallback, useContext, useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { DoesNotExist, useApi } from '../api'
import { useErrors } from '../errors'
import { useHistory } from 'react-router-dom'
import { getUrl } from '../nav/Routes'

export const uploadContext = React.createContext()

export function useUploadContext() {
  return useContext(uploadContext)
}

const UploadContext = React.memo(function UploadContext({uploadId, children}) {
  const {api, user} = useApi()
  const {raiseError} = useErrors()
  const history = useHistory()

  const [pagination, setPagination] = useState({
    page_size: 5, page: 1, order: 'asc', order_by: 'process_status'
  })
  const [deleteClicked, setDeleteClicked] = useState(false)
  const [data, setData] = useState(null)
  const [apiData, setApiData] = useState(null)
  const [error, setError] = useState(null)
  const upload = data?.upload
  const hasUpload = !!upload
  const setUpload = useMemo(() => (upload) => {
    setData(data => ({...data, upload: upload}))
  }, [setData])
  const isProcessing = upload?.process_running

  const fetchData = useCallback(() => () => {
    api.get(`/uploads/${uploadId}/entries`, pagination, {returnRequest: true})
      .then(apiData => {
        setApiData(apiData)
        setData(apiData.response)
      })
      .catch((error) => {
        if (error instanceof DoesNotExist && deleteClicked) {
          history.push(getUrl('uploads', new URL(window.location.href).pathname))
          return
        }
        if (error.apiMessage === 'Page out of range requested.') {
          // Can happen if entries have been deleted and the page is no longer in range
          pagination.page = 1
          setPagination(pagination)
        } else if (!hasUpload && error.apiMessage) {
          setError(error.apiMessage)
        } else {
          raiseError(error)
        }
      })
  }, [api, hasUpload, uploadId, pagination, setPagination, deleteClicked, raiseError, setData, setApiData, history])

  // constant fetching of upload data when necessary
  useEffect(() => {
    if (isProcessing) {
      const interval = setInterval(fetchData(), 1000)
      return () => clearInterval(interval)
    }
  }, [fetchData, isProcessing])

  // initial fetching of upload data
  useEffect(fetchData(), [fetchData])

  const viewers = upload?.viewers
  const writers = upload?.writers
  const isViewer = user && viewers?.includes(user.sub)
  const isWriter = user && writers?.includes(user.sub)

  const contextValue = useMemo(() => ({
    uploadId: uploadId,
    upload: upload,
    hasUpload: hasUpload,
    setUpload: setUpload,
    setPagination: setPagination,
    pagination: pagination,
    data: data,
    apiData: apiData,
    isViewer: isViewer,
    isWriter: isWriter,
    update: fetchData,
    error: error,
    deleteClicked: deleteClicked,
    setDeleteClicked: setDeleteClicked,
    isProcessing: isProcessing
  }), [
    uploadId, upload, hasUpload, setUpload, setPagination, pagination, data, apiData,
    isViewer, isWriter, fetchData, error, deleteClicked, setDeleteClicked, isProcessing
  ])

  return <uploadContext.Provider value={contextValue}>
    {children}
  </uploadContext.Provider>
})
UploadContext.propTypes = {
  uploadId: PropTypes.string.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default UploadContext
