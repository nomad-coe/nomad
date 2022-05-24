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
import React, { useContext, useRef } from 'react'
import PropTypes from 'prop-types'
import { useApi, DoesNotExist } from './api'
import { useErrors } from './errors'

function addSubscription(storeObj, cb, options) {
  storeObj.subscriptions.push({cb, ...options})
}

function removeSubscription(store, key, cb) {
  const storeObj = store[key]
  if (storeObj) {
    storeObj.subscriptions = storeObj.subscriptions.filter(subscription => subscription.cb !== cb)
  }
}

export const dataStoreContext = React.createContext()

export function useDataStoreContext() {
  return useContext(dataStoreContext)
}

const DataStore = React.memo(({children}) => {
  const {raiseError} = useErrors()
  const {api, user} = useApi()
  const userRef = useRef()
  userRef.current = user // The logged in user may change during component lifecycle

  const uploadStore = useRef({}) // The upload store objects

  function getUpload(uploadId) {
    if (!uploadId) return undefined
    let uploadStoreObj = uploadStore.current[uploadId]
    if (!uploadStoreObj) {
      // Creates an initial, empty upload store object.
      uploadStoreObj = {
        uploadId: uploadId, // ReadOnly
        deletionRequested: false, // Writable - If this upload has been sent for deletion
        upload: undefined, // Writeable - The last upload proc data fetched.
        entries: undefined, // ReadOnly - The last list of entries fetched by the store (when subscribing to an entry page).
        apiData: undefined, // ReadOnly - Object with the last api request and response fetched by the store (when subscribing to an entry page).
        pagination: {
          page_size: 5, page: 1, order: 'asc', order_by: 'process_status'
        }, // Writeable - Value for pagination to use when an entries page is requested

        // ReadOnly - Derived values (computed when calling updateUpload)
        hasUpload: false,
        isProcessing: false,
        isViewer: false,
        isWriter: false,
        isEditable: false,

        // ReadOnly - Managed by the store
        error: undefined, // If we had an api error from the last store refresh
        subscriptions: [],
        refreshing: false
      }
      uploadStore.current[uploadId] = uploadStoreObj
    }
    return uploadStoreObj
  }

  function subscribeToUpload(uploadId, cb, requireUpload, requireEntriesPage) {
    // Subscribes the callback cb to an upload, and returns a function to be called to unsubscribe.
    // Designed to be used in useEffect. The callback will be called when the stored value changes.
    if (!uploadId) return undefined
    if (requireUpload === undefined || requireEntriesPage === undefined) {
      throw Error('Store error: missing upload subscription parameter')
    }
    const uploadStoreObj = getUpload(uploadId)
    addSubscription(uploadStoreObj, cb, {requireUpload, requireEntriesPage})
    initiateUploadRefreshIfNeeded(uploadId)
    return function unsubscriber() { removeSubscription(uploadStore.current, uploadId, cb) }
  }

  function updateUpload(uploadId, dataToUpdate) {
    // Updates the specified data and notifies all subscribers
    const oldStoreObj = getUpload(uploadId)
    const newStoreObj = {...oldStoreObj, ...dataToUpdate}
    // Compute derived values
    const user = userRef.current
    newStoreObj.hasUpload = !!newStoreObj.upload
    newStoreObj.isProcessing = !!newStoreObj.upload?.process_running
    if (dataToUpdate.upload || dataToUpdate.entries) {
      dataToUpdate.error = undefined // Updating upload or entries -> reset error
    }
    const viewers = newStoreObj.upload?.viewers
    const writers = newStoreObj.upload?.writers
    newStoreObj.isViewer = user && viewers?.includes(user.sub)
    newStoreObj.isWriter = user && writers?.includes(user.sub)
    newStoreObj.isEditable = newStoreObj.isWriter && !newStoreObj.upload.published

    // Update the store
    uploadStore.current[uploadId] = newStoreObj

    // Notify subscribers
    for (const subscription of newStoreObj.subscriptions) {
      try {
        subscription.cb(oldStoreObj, newStoreObj)
      } catch (error) {
        console.error('DataStore: failed to notify subscriber: ' + error.message)
      }
    }
    // Report any unexpected api errors to the user
    if (newStoreObj.error && (!oldStoreObj.error || oldStoreObj.error?.apiMessage !== newStoreObj.error?.apiMessage)) {
      if (newStoreObj.error instanceof DoesNotExist && newStoreObj.deletionRequested) {
        // Expected to happen when user deletes an upload, ignore
      } else {
        // Unexpected error occured when refreshing
        raiseError(newStoreObj.error)
      }
    }
    // Possibly, start a refresh job
    initiateUploadRefreshIfNeeded(uploadId)
  }

  function uploadDataRequired(uploadStoreObj) {
    // Determine, by looking at the subscriptions, which data the store needs to fetch for an upload
    let requireUpload = false
    let requireEntriesPage = false
    for (const subscription of uploadStoreObj.subscriptions) {
      requireUpload |= subscription.requireUpload
      requireEntriesPage |= subscription.requireEntriesPage
    }
    return [requireUpload, requireEntriesPage]
  }

  async function refreshUpload(uploadId) {
    // Used internally by the store to refresh an upload store obj with data from the API.
    const uploadStoreObj = getUpload(uploadId)
    const [requireUpload, requireEntriesPage] = uploadDataRequired(uploadStoreObj)
    if (!requireUpload && !requireEntriesPage) return
    // There's something that can be refreshed
    uploadStoreObj.refreshing = true
    const currentPagination = uploadStoreObj.pagination

    const apiCall = requireEntriesPage
      ? api.get(`/uploads/${uploadId}/entries`, currentPagination, {returnRequest: true})
      : api.get(`/uploads/${uploadId}`)

    apiCall.then(apiData => {
      const dataToUpdate = requireEntriesPage
        ? {error: undefined, refreshing: false, upload: apiData.response?.upload, entries: apiData.response?.data, apiData: apiData, pagination: currentPagination}
        : {error: undefined, refreshing: false, upload: apiData.data, entries: undefined, apiData: undefined}
      updateUpload(uploadId, dataToUpdate)
    }).catch((error) => {
      if (requireEntriesPage && error.apiMessage === 'Page out of range requested.') {
        // Special case: can happen if entries have been deleted and the page we were on is no longer in range
        if (currentPagination && currentPagination.page !== 1) {
          // Rather than sending an update to all subscribers with an error, we first try
          // jumping to page 1 (will probably solve the problem)
          getUpload(uploadId).pagination.page = 1
          refreshUpload(uploadId)
          return
        }
      }
      updateUpload(uploadId, {error: error, refreshing: false})
    })
  }

  async function requestRefreshUpload(uploadId) {
    // Used from child components or asyncronously to request a refresh of the upload store object
    const uploadStoreObj = getUpload(uploadId)
    if (!uploadStoreObj.refreshing) {
      // Refresh is not already in progress
      refreshUpload(uploadId)
    }
  }

  async function initiateUploadRefreshIfNeeded(uploadId) {
    // Called when an upload is first subscribed to, and whenever it is updated using updateUpload
    let uploadStoreObj = getUpload(uploadId)
    if (uploadStoreObj.refreshing) return // refresh already in progress
    if (uploadStoreObj.isProcessing) {
      // Upload is processing
      uploadStoreObj.refreshing = true // Signal start of a refresh
      await new Promise(resolve => setTimeout(resolve, 1000)) // wait one sec
      uploadStoreObj = getUpload(uploadId)
    }
    // Determine if a refresh is needed or not
    const [requireUpload, requireEntriesPage] = uploadDataRequired(uploadStoreObj)
    const uploadDataMissing = requireUpload && !uploadStoreObj.upload
    const entryDataMissing = requireEntriesPage && !uploadStoreObj.entries
    const pag = uploadStoreObj.pagination
    const pagIs = uploadStoreObj.apiData?.response?.pagination
    const wrongPagination = requireEntriesPage && (pagIs?.page !== pag?.page || pagIs?.page_size !== pag.page_size)
    if (!uploadStoreObj.error && (uploadDataMissing || entryDataMissing || wrongPagination || uploadStoreObj.isProcessing)) {
      // Need to fetch data from the api
      refreshUpload(uploadId)
    } else {
      uploadStoreObj.refreshing = false
    }
  }

  const contextValue = {
    getUpload: getUpload,
    subscribeToUpload: subscribeToUpload,
    updateUpload: updateUpload,
    requestRefreshUpload: requestRefreshUpload
  }

  return <dataStoreContext.Provider value={contextValue}>
    {children}
  </dataStoreContext.Provider>
})
DataStore.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}
export default DataStore
