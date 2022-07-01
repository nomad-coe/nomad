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

export function useDataStore() {
  return useContext(dataStoreContext)
}

const DataStore = React.memo(({children}) => {
  const {raiseError} = useErrors()
  const {api, user} = useApi()
  const userRef = useRef()
  userRef.current = user // The logged in user may change during component lifecycle

  const uploadStore = useRef({}) // The upload store objects
  const entryStore = useRef({}) // The entry store objects

  /**
   * Gets an upload object from the store, creating it if it doesn't exist (in which case
   * an object with default values, mostly undefined or nulls, will be returned). Note, it
   * does not cause the store to fetch any data; for that, a subscription is necessary.
   */
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
        isRefreshing: false,
        refreshOptions: null // The options used at the last store refresh
      }
      uploadStore.current[uploadId] = uploadStoreObj
    }
    return uploadStoreObj
  }

  /**
   * Gets an upload from the store asychronously, waiting for the store to refresh if needed.
   * If the store needs to be refreshed, the call will set up a temporary subscription, which
   * will be removed when finished. Note that to keep the object in the store and watch it
   * also after this call has ended, you need to set up a subscription.
   */
  async function getUploadAsync(uploadId, requireUpload, requireEntriesPage) {
    if (!uploadId) return undefined
    const uploadStoreObj = getUpload(uploadId)
    if (uploadRefreshSatisfiesOptions(uploadStoreObj, requireUpload, requireEntriesPage)) {
      return uploadStoreObj // Store has already been refreshed with the required options
    }
    // Set up temporary subscription and wait for it to refresh
    return new Promise(function(resolve, reject) {
      function cb(oldStoreObj, newStoreObj) {
        if (uploadRefreshSatisfiesOptions(newStoreObj, requireUpload, requireEntriesPage)) {
          removeSubscription(uploadStore.current, uploadId, cb)
          resolve(newStoreObj)
        }
      }
      subscribeToUpload(uploadId, cb, requireUpload, requireEntriesPage)
    })
  }

  /**
   * Subscribes the callback cb to an upload, and returns a function to be called to unsubscribe.
   * Typically used in useEffect. The callback will be called when the store value changes.
   */
  function subscribeToUpload(uploadId, cb, requireUpload, requireEntriesPage) {
    if (!uploadId) return undefined
    if (requireUpload === undefined || requireEntriesPage === undefined) {
      throw Error('Store error: missing upload subscription parameter')
    }
    const uploadStoreObj = getUpload(uploadId)
    addSubscription(uploadStoreObj, cb, {requireUpload, requireEntriesPage})
    initiateUploadRefreshIfNeeded(uploadId)
    return function unsubscriber() { removeSubscription(uploadStore.current, uploadId, cb) }
  }

  /**
   * Updates the store upload with the specified data and notifies all subscribers.
   */
  function updateUpload(uploadId, dataToUpdate) {
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
    for (const subscription of [...newStoreObj.subscriptions]) {
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

  function uploadOptions(uploadStoreObj) {
    // Internal use: Determine the options of the upload, based on the subscriptions
    let requireUpload = false
    let requireEntriesPage = false
    for (const subscription of uploadStoreObj.subscriptions) {
      requireUpload = requireUpload || subscription.requireUpload
      requireEntriesPage = requireEntriesPage || subscription.requireEntriesPage
    }
    return {requireUpload, requireEntriesPage}
  }

  function uploadRefreshSatisfiesOptions(uploadStoreObj, requireUpload, requireEntriesPage) {
    // Internal use: Determine if the last store refresh satisfies the given options.
    const refreshOptions = uploadStoreObj?.refreshOptions
    if (refreshOptions) {
      return (refreshOptions.requireUpload || !requireUpload) || (refreshOptions.requireEntriesPage || !requireEntriesPage)
    }
    return false
  }

  async function refreshUpload(uploadId) {
    // Internal use: refresh an upload store obj with data from the API.
    const uploadStoreObj = getUpload(uploadId)
    const refreshOptions = uploadOptions(uploadStoreObj)
    const {requireUpload, requireEntriesPage} = refreshOptions
    if (!requireUpload && !requireEntriesPage) return
    // There's something that can be refreshed
    uploadStoreObj.isRefreshing = true
    const currentPagination = uploadStoreObj.pagination

    const apiCall = requireEntriesPage
      ? api.get(`/uploads/${uploadId}/entries`, currentPagination, {returnRequest: true})
      : api.get(`/uploads/${uploadId}`)

    apiCall.then(apiData => {
      const dataToUpdate = requireEntriesPage
        ? {error: undefined, isRefreshing: false, upload: apiData.response?.upload, entries: apiData.response?.data, apiData, pagination: currentPagination, refreshOptions}
        : {error: undefined, isRefreshing: false, upload: apiData.data, entries: undefined, apiData: undefined, refreshOptions}
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
      updateUpload(uploadId, {error: error, isRefreshing: false, refreshOptions})
    })
  }

  /**
   * Use to nicely request a refresh of the upload store object.
   */
  function requestRefreshUpload(uploadId) {
    const uploadStoreObj = getUpload(uploadId)
    if (!uploadStoreObj.isRefreshing) {
      // Refresh is not already in progress
      refreshUpload(uploadId)
    }
  }

  async function initiateUploadRefreshIfNeeded(uploadId) {
    // Internal use: check if a refresh of the store is needed, and if so, initiate it.
    let uploadStoreObj = getUpload(uploadId)
    if (uploadStoreObj.isRefreshing) return // refresh already in progress
    if (uploadStoreObj.isProcessing) {
      // Upload is processing
      uploadStoreObj.isRefreshing = true // Signal start of a refresh
      await new Promise(resolve => setTimeout(resolve, 1000)) // wait one sec
      uploadStoreObj = getUpload(uploadId)
    }
    // Determine if a refresh is needed or not
    const {requireUpload, requireEntriesPage} = uploadOptions(uploadStoreObj)
    const uploadDataMissing = requireUpload && !uploadStoreObj.upload
    const entryDataMissing = requireEntriesPage && !uploadStoreObj.entries
    const pag = uploadStoreObj.pagination
    const pagIs = uploadStoreObj.apiData?.response?.pagination
    const wrongPagination = requireEntriesPage && (pagIs?.page !== pag?.page || pagIs?.page_size !== pag.page_size)
    if (!uploadStoreObj.error && (uploadDataMissing || entryDataMissing || wrongPagination || uploadStoreObj.isProcessing)) {
      // Need to fetch data from the api
      refreshUpload(uploadId)
    } else {
      uploadStoreObj.isRefreshing = false
    }
  }

  /**
   * Gets an entry object from the store, creating it if it doesn't exist (in which case
   * an object with default values, mostly undefined or nulls, will be returned). Note, it
   * does not cause the store to fetch any data; for that, a subscription is necessary.
   */
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
        isRefreshing: false,
        refreshOptions: null, // The options used at the last store refresh

        // Convenience methods
        handleArchiveChanged: () => { handleArchiveChanged(entryId) },
        saveArchive: () => { return saveArchive(entryId) },
        reload: () => { requestRefreshEntry(entryId) }
      }
      entryStore.current[entryId] = entryStoreObj
    }
    return entryStoreObj
  }

  /**
   * Gets an entry from the store asychronously, waiting for the store to refresh if needed.
   * If the store needs to be refreshed, the call will set up a temporary subscription, which
   * will be removed when finished. Note that to keep the object in the store and watch it
   * also after this call has ended, you need to set up a subscription.
   */
  async function getEntryAsync(entryId, requireMetadata, requireArchive) {
    if (!entryId) return undefined
    const entryStoreObj = getEntry(entryId)
    if (entryRefreshSatisfiesOptions(entryStoreObj, requireMetadata, requireArchive)) {
      return entryStoreObj // Store has already been refreshed with the required options
    }
    // Set up temporary subscription and wait for it to refresh
    return new Promise(function(resolve, reject) {
      function cb(oldStoreObj, newStoreObj) {
        if (entryRefreshSatisfiesOptions(newStoreObj, requireMetadata, requireArchive)) {
          removeSubscription(entryStore.current, entryId, cb)
          resolve(newStoreObj)
        }
      }
      subscribeToEntry(entryId, cb, requireMetadata, requireArchive, false)
    })
  }

  /**
   * Subscribes the callback cb to an entry, and returns a function to be called to unsubscribe.
   * Typically used in useEffect. The callback will be called when the store value changes.
   */
  function subscribeToEntry(entryId, cb, requireMetadata, requireArchive, editableIfPossible) {
    if (!entryId) return undefined
    if (requireMetadata === undefined || requireArchive === undefined || editableIfPossible === undefined) {
      throw Error('Store error: missing entry subscription parameter')
    }
    const entryStoreObj = getEntry(entryId)
    addSubscription(entryStoreObj, cb, {requireMetadata, requireArchive, editableIfPossible})
    initiateEntryRefreshIfNeeded(entryId)
    return function unsubscriber() { removeSubscription(entryStore.current, entryId, cb) }
  }

  /**
   * Updates the store entry with the specified data and notifies all subscribers.
   */
  function updateEntry(entryId, dataToUpdate) {
    const oldStoreObj = getEntry(entryId)
    const newStoreObj = {...oldStoreObj, ...dataToUpdate}
    // Compute derived values not set by the refreshEntry method
    newStoreObj.exists = newStoreObj?.error?.name !== 'DoesNotExist'
    newStoreObj.archiveHasChanges = newStoreObj.archiveVersion !== newStoreObj.savedArchiveVersion

    // Update the store
    entryStore.current[entryId] = newStoreObj

    // Notify subscribers
    for (const subscription of [...newStoreObj.subscriptions]) {
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

  function entryOptions(entryStoreObj) {
    // Internal use: Determine the options of the entry, based on the subscriptions
    let requireMetadata = false
    let requireArchive = false
    let editableIfPossible = false
    for (const subscription of entryStoreObj.subscriptions) {
      requireMetadata = requireMetadata || subscription.requireMetadata
      requireArchive = requireArchive || subscription.requireArchive
      editableIfPossible = editableIfPossible || subscription.editableIfPossible
    }
    return {requireMetadata, requireArchive, editableIfPossible}
  }

  function entryRefreshSatisfiesOptions(entryStoreObj, requireMetadata, requireArchive) {
    // Internal use: Determine if the last store refresh satisfies the given options.
    const refreshOptions = entryStoreObj?.refreshOptions
    if (refreshOptions) {
      return (refreshOptions.requireMetadata || !requireMetadata) || (refreshOptions.requireArchive || !requireArchive)
    }
    return false
  }

  async function refreshEntry(entryId) {
    // Internal use: refresh an entry store obj with data from the API.
    const entryStoreObj = getEntry(entryId)
    const {requireMetadata, requireArchive, editableIfPossible} = entryOptions(entryStoreObj)
    if (!requireMetadata && !requireArchive) return
    // There's something that can be refreshed
    entryStoreObj.isRefreshing = true
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
    dataToUpdate.isRefreshing = false
    updateEntry(entryId, dataToUpdate)
  }

  /**
   * Use to nicely request a refresh of the entry store object.
   */
  function requestRefreshEntry(entryId) {
    const entryStoreObj = getEntry(entryId)
    if (!entryStoreObj.isRefreshing) {
      // Refresh is not already in progress
      refreshEntry(entryId)
    }
  }

  async function initiateEntryRefreshIfNeeded(entryId) {
    // Internal use: check if a refresh of the store is needed, and if so, initiate it.
    let entryStoreObj = getEntry(entryId)
    if (entryStoreObj.isRefreshing) return // refresh already in progress
    if (entryStoreObj.isProcessing) {
      // Entry is processing
      entryStoreObj.isRefreshing = true // Signal start of a refresh
      await new Promise(resolve => setTimeout(resolve, 1000)) // wait one sec
      entryStoreObj = getEntry(entryId)
    }
    // Determine if a refresh is needed or not
    const {requireMetadata, requireArchive} = entryOptions(entryStoreObj)
    const metadataMissing = requireMetadata && !entryStoreObj.metadata
    const archiveMissing = requireArchive && !entryStoreObj.archive
    if (!entryStoreObj.error && (metadataMissing || archiveMissing || entryStoreObj.isProcessing)) {
      // Need to fetch data from the api
      refreshEntry(entryId)
    } else {
      entryStoreObj.isRefreshing = false
    }
  }

  /**
   * Used to save the archive and trigger a store refresh.
   */
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
      return new Promise((resolve, reject) => {
        api.put(`/uploads/${uploadId}/raw/${path}?file_name=${fileName}&wait_for_processing=true&entry_hash=${archive.metadata.entry_hash}`, newArchive)
          .then(response => {
            requestRefreshEntry(entryId)
          })
          .catch(error => {
            if (error?.status === 409) {
              reject(error)
            } else {
              raiseError(error)
            }
          })
        updateEntry(entryId, {savedArchiveVersion: archiveVersion})
      })
    }
  }

  /**
   * Call to signal that the archive has been manually edited.
   */
  function handleArchiveChanged(entryId) {
    const {archiveVersion} = getEntry(entryId)
    updateEntry(entryId, {archiveVersion: archiveVersion + 1})
  }

  const contextValue = {
    getUpload,
    getUploadAsync,
    subscribeToUpload,
    updateUpload,
    requestRefreshUpload,
    getEntry,
    getEntryAsync,
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
