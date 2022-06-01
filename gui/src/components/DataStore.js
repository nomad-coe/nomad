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
    if (!storeObj.subscriptions.length) {
      // No subscribers. Wait for 5 seconds, if there are still no subscribers, clear the cache
      setTimeout(() => {
        if (store) {
          const storeObj = store[key]
          if (storeObj && !storeObj.subscriptions.length) {
            delete store[key]
          }
        }
      }, 5000)
    }
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
  const entryStore = useRef({}) // The entry store objects

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

  function uploadSubscriptionOptions(uploadStoreObj) {
    // Determine the subscription options of the upload
    let requireUpload = false
    let requireEntriesPage = false
    for (const subscription of uploadStoreObj.subscriptions) {
      requireUpload |= subscription.requireUpload
      requireEntriesPage |= subscription.requireEntriesPage
    }
    return [!!requireUpload, !!requireEntriesPage]
  }

  async function refreshUpload(uploadId) {
    // Used internally by the store to refresh an upload store obj with data from the API.
    const uploadStoreObj = getUpload(uploadId)
    const [requireUpload, requireEntriesPage] = uploadSubscriptionOptions(uploadStoreObj)
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
    const [requireUpload, requireEntriesPage] = uploadSubscriptionOptions(uploadStoreObj)
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

  function getEntry(entryId) {
    if (!entryId) return undefined
    let entryStoreObj = entryStore.current[entryId]
    if (!entryStoreObj) {
      // Creates an initial, empty entry store object.
      entryStoreObj = {
        entryId: entryId, // ReadOnly
        uploadId: undefined, // ReadOnly - fetched by the store
        metadata: undefined, // ReadOnly - fetched by the store
        metadataApiData: undefined, // ReadOnly - fetched by the store
        archive: undefined, // Modifiable object - fetched by the store, but the object *content* can be changed when editing.
        archiveApiData: undefined, // Modifiable object - fetched by the store, but the "archive" key is the same object as above, i.e. modifiable.

        // ReadOnly - Derived or managed by the store
        exists: true,
        isProcessing: false,
        editable: false,
        archiveVersion: 0, // Used to keep track of manual edits
        savedArchiveVersion: 0,
        archiveHasChanges: false,

        error: undefined, // If we had an api error from the last store refresh
        subscriptions: [],
        refreshing: false,

        // Convenience methods
        handleArchiveChanged: () => { handleArchiveChanged(entryId) },
        saveArchive: () => { saveArchive(entryId) }
      }
      entryStore.current[entryId] = entryStoreObj
    }
    return entryStoreObj
  }

  function subscribeToEntry(entryId, cb, requireMetadata, requireArchive, editableIfPossible) {
    // Subscribes the callback cb to an entry, and returns a function to be called to unsubscribe.
    // Designed to be used in useEffect. The callback will be called when the stored value changes.
    if (!entryId) return undefined
    if (requireMetadata === undefined || requireArchive === undefined || editableIfPossible === undefined) {
      throw Error('Store error: missing entry subscription parameter')
    }
    const entryStoreObj = getEntry(entryId)
    addSubscription(entryStoreObj, cb, {requireMetadata, requireArchive, editableIfPossible})
    initiateEntryRefreshIfNeeded(entryId)
    return function unsubscriber() { removeSubscription(entryStore.current, entryId, cb) }
  }

  function updateEntry(entryId, dataToUpdate) {
    // Updates the specified data and notifies all subscribers
    const oldStoreObj = getEntry(entryId)
    const newStoreObj = {...oldStoreObj, ...dataToUpdate}
    // Compute derived values not set by the refreshEntry method
    newStoreObj.exists = newStoreObj?.error?.name !== 'DoesNotExist'
    newStoreObj.archiveHasChanges = newStoreObj.archiveVersion !== newStoreObj.savedArchiveVersion

    // Update the store
    entryStore.current[entryId] = newStoreObj

    // Notify subscribers
    for (const subscription of newStoreObj.subscriptions) {
      try {
        subscription.cb(oldStoreObj, newStoreObj)
      } catch (error) {
        console.error('DataStore: failed to notify subscriber: ' + error.message)
      }
    }
    // Report any unexpected api errors to the user
    if (newStoreObj.exists && newStoreObj.error && (!oldStoreObj.error || oldStoreObj.error?.apiMessage !== newStoreObj.error?.apiMessage)) {
      // Unexpected error occured when refreshing
      raiseError(newStoreObj.error)
    }
    // Possibly, start a refresh job
    initiateEntryRefreshIfNeeded(entryId)
  }

  function entrySubscriptionOptions(entryStoreObj) {
    // Determine the subscription options of the entry
    let requireMetadata = false
    let requireArchive = false
    let editableIfPossible = false
    for (const subscription of entryStoreObj.subscriptions) {
      requireMetadata |= subscription.requireMetadata
      requireArchive |= subscription.requireArchive
      editableIfPossible |= subscription.editableIfPossible
    }
    return [!!requireMetadata, !!requireArchive, !!editableIfPossible]
  }

  async function refreshEntry(entryId) {
    // Used internally by the store to refresh an entry store obj with data from the API.
    const entryStoreObj = getEntry(entryId)
    const [requireMetadata, requireArchive, editableIfPossible] = entrySubscriptionOptions(entryStoreObj)
    if (!requireMetadata && !requireArchive) return
    // There's something that can be refreshed
    entryStoreObj.refreshing = true
    let dataToUpdate = {}
    try {
      const metadataApiData = await api.get(`/entries/${entryId}`, null, {returnRequest: true})
      const metadata = metadataApiData?.response?.data
      const uploadId = metadata?.upload_id
      const user = userRef.current
      const isWriter = user && metadata?.writers && metadata.writers.find(u => u.user_id === user.sub)
      const isEditableArchive = metadata && !metadata.published && metadata.quantities && metadata.quantities.includes('data')
      const editable = editableIfPossible && isWriter && isEditableArchive
      const isProcessing = !!metadata?.process_running
      dataToUpdate = {metadataApiData, metadata, uploadId, editable, isProcessing, error: undefined}
      if (requireArchive) {
        const required = editable
          ? '*'
          : {
            'resolve-inplace': false,
            metadata: '*',
            data: '*',
            definitions: '*',
            results: {
              material: '*',
              method: '*',
              properties: {
                structures: '*',
                electronic: 'include-resolved',
                mechanical: 'include-resolved',
                spectroscopy: 'include-resolved',
                vibrational: 'include-resolved',
                thermodynamic: 'include-resolved',
                // For geometry optimizations we require only the energies.
                // Trajectory, optimized structure, etc. are unnecessary.
                geometry_optimization: {
                  energies: 'include-resolved'
                }
              }
            }
          }
        const archiveApiData = await api.post(
          `/entries/${entryId}/archive/query`, {required}, {returnRequest: true, jsonResponse: true})
        const archive = archiveApiData?.response?.data?.archive
        if (archive) {
          archive.processing_logs = undefined
        }
        dataToUpdate.archiveApiData = archiveApiData
        dataToUpdate.archive = archive
      }
    } catch (error) {
      dataToUpdate.error = error
    }
    dataToUpdate.refreshing = false
    updateEntry(entryId, dataToUpdate)
  }

  async function requestRefreshEntry(entryId) {
    // Used from child components or asyncronously to request a refresh of the entry store object
    const entryStoreObj = getEntry(entryId)
    if (!entryStoreObj.refreshing) {
      // Refresh is not already in progress
      refreshEntry(entryId)
    }
  }

  async function initiateEntryRefreshIfNeeded(entryId) {
    // Called when an entry is first subscribed to, and whenever it is updated using updateEntry
    let entryStoreObj = getEntry(entryId)
    if (entryStoreObj.refreshing) return // refresh already in progress
    if (entryStoreObj.isProcessing) {
      // Entry is processing
      entryStoreObj.refreshing = true // Signal start of a refresh
      await new Promise(resolve => setTimeout(resolve, 1000)) // wait one sec
      entryStoreObj = getEntry(entryId)
    }
    // Determine if a refresh is needed or not
    const [requireMetadata, requireArchive] = entrySubscriptionOptions(entryStoreObj)
    const metadataMissing = requireMetadata && !entryStoreObj.metadata
    const archiveMissing = requireArchive && !entryStoreObj.archive
    if (!entryStoreObj.error && (metadataMissing || archiveMissing || entryStoreObj.isProcessing)) {
      // Need to fetch data from the api
      refreshEntry(entryId)
    } else {
      entryStoreObj.refreshing = false
    }
  }

  function saveArchive(entryId) {
    const {uploadId, metadata, archive, archiveVersion} = getEntry(entryId)
    const {mainfile} = metadata
    if (uploadId) {
      const separatorIndex = mainfile.lastIndexOf('/')
      const path = mainfile.substring(0, separatorIndex + 1)
      const fileName = mainfile.substring(separatorIndex + 1)
      const newArchive = {...archive}
      delete newArchive.metadata
      delete newArchive.results
      delete newArchive.processing_logs
      api.put(`/uploads/${uploadId}/raw/${path}?file_name=${fileName}&wait_for_processing=true`, newArchive)
        .then(response => {
          requestRefreshEntry(entryId)
        })
        .catch(raiseError)
      updateEntry(entryId, {savedArchiveVersion: archiveVersion})
    }
  }

  function handleArchiveChanged(entryId) {
    const {archiveVersion} = getEntry(entryId)
    updateEntry(entryId, {archiveVersion: archiveVersion + 1})
  }

  const contextValue = {
    getUpload,
    subscribeToUpload,
    updateUpload,
    requestRefreshUpload,
    getEntry,
    subscribeToEntry
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
