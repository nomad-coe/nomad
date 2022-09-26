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
import React, { useContext, useRef, useState, useCallback, useEffect } from 'react'
import PropTypes from 'prop-types'
import { useApi, DoesNotExist } from './api'
import { useErrors } from './errors'
import { apiBase } from '../config'
import { refType, parseNomadUrl, createEntryUrl, systemMetainfoUrl } from '../utils'
import { getUrlFromDefinition, Metainfo } from './archive/metainfo'
import currentSystemMetainfoData from '../metainfo'
import YAML from 'yaml'

function addSubscription(storeObj, cb, options) {
  storeObj._subscriptions.push({cb, ...options})
  storeObj._subscribeCounter = (storeObj._subscribeCounter || 0) + 1
}

function removeSubscription(store, key, cb) {
  const storeObj = store[key]
  if (storeObj) {
    storeObj._subscriptions = storeObj._subscriptions.filter(subscription => subscription.cb !== cb)
    if (!storeObj._subscriptions.length) {
      // No subscribers. Wait for 5 seconds, then check again
      const oldSubscribeCount = storeObj._subscribeCounter
      setTimeout(() => {
        if (store) {
          const storeObj = store[key]
          if (storeObj && storeObj._subscribeCounter === oldSubscribeCount) {
            // Noone has subscribed to it during this time window. Delete it.
            delete store[key]
          }
        }
      }, 5000)
    }
  }
}

const metainfoArchiveFilter = Object.freeze({
  // Archive data required to create a Metainfo object
  definitions: '*',
  metadata: '*'
})

export const dataStoreContext = React.createContext()

export function useDataStore() {
  return useContext(dataStoreContext)
}

const DataStore = React.memo(({children}) => {
  const {raiseError} = useErrors()
  const {api, user} = useApi()
  const userRef = useRef()
  userRef.current = user // The logged in user may change during component lifecycle
  const selectedEntry = useRef() // stores identity of the selected entry (used to determine if editable)

  const uploadStore = useRef({}) // The upload store objects
  const entryStore = useRef({}) // The entry store objects
  const metainfoDataStore = useRef({}) // The metainfo data store objects
  const externalInheritanceCache = useRef({}) // Used to keep track of inheritance between metainfos

  /**
   * Gets an upload object from the store, creating it if it doesn't exist (in which case
   * an object with default values, mostly undefined or nulls, will be returned). Note, it
   * does not cause the store to fetch any data; for that, a subscription is necessary.
   */
  function getUpload(installationUrl, uploadId) {
    if (!uploadId) return undefined
    if (installationUrl !== apiBase) throw new Error('Fetching uploads from external installations is not yet supported')
    let uploadStoreObj = uploadStore.current[uploadId]
    if (!uploadStoreObj) {
      // Creates an initial, empty upload store object.
      uploadStoreObj = {
        installationUrl, // ReadOnly
        uploadId, // ReadOnly
        isExternal: installationUrl !== apiBase, // ReadOnly
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
        isRefreshing: false,
        refreshOptions: null, // The options used at the last store refresh
        _subscriptions: [],

        // Convenience methods
        updateUpload: (dataToUpdate) => { updateUpload(installationUrl, uploadId, dataToUpdate) },
        requestRefreshUpload: () => { requestRefreshUpload(installationUrl, uploadId) }
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
  async function getUploadAsync(installationUrl, uploadId, requireUpload, requireEntriesPage) {
    if (!uploadId) return undefined
    const uploadStoreObj = getUpload(installationUrl, uploadId)
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
      subscribeToUpload(installationUrl, uploadId, cb, requireUpload, requireEntriesPage)
    })
  }

  /**
   * Subscribes the callback cb to an upload, and returns a function to be called to unsubscribe.
   * Typically used in useEffect. The callback will be called when the store value changes.
   */
  function subscribeToUpload(installationUrl, uploadId, cb, requireUpload, requireEntriesPage) {
    if (!uploadId) return undefined
    if (requireUpload === undefined || requireEntriesPage === undefined) {
      throw Error('Store error: missing upload subscription parameter')
    }
    const uploadStoreObj = getUpload(installationUrl, uploadId)
    addSubscription(uploadStoreObj, cb, {requireUpload, requireEntriesPage})
    initiateUploadRefreshIfNeeded(installationUrl, uploadId)
    return function unsubscriber() { removeSubscription(uploadStore.current, uploadId, cb) }
  }

  /**
   * Updates the store upload with the specified data and notifies all subscribers.
   */
  function updateUpload(installationUrl, uploadId, dataToUpdate) {
    if (installationUrl !== apiBase) throw new Error('Cannot update external upload')
    const oldStoreObj = getUpload(installationUrl, uploadId)
    const newStoreObj = {...oldStoreObj, ...dataToUpdate}
    // Compute derived values
    const user = userRef.current
    newStoreObj.hasUpload = !!newStoreObj.upload
    newStoreObj.isProcessing = !!newStoreObj.upload?.process_running
    if (dataToUpdate.upload || dataToUpdate.entries) {
      dataToUpdate.error = undefined // Updating upload or entries -> reset error
    }
    if (dataToUpdate?.upload?.current_process === 'delete_upload') {
      newStoreObj.deletionRequested = true // Will treat subsequent 404 errors
    }
    const viewers = newStoreObj.upload?.viewers
    const writers = newStoreObj.upload?.writers
    newStoreObj.isViewer = user && viewers?.includes(user.sub)
    newStoreObj.isWriter = user && writers?.includes(user.sub)
    newStoreObj.isEditable = newStoreObj.isWriter && !newStoreObj.upload.published

    // Update the store
    uploadStore.current[uploadId] = newStoreObj

    // Notify subscribers
    for (const subscription of [...newStoreObj._subscriptions]) {
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
    initiateUploadRefreshIfNeeded(installationUrl, uploadId)
  }

  async function refreshUpload(installationUrl, uploadId) {
    // Internal use: refresh an upload store obj with data from the API.
    const uploadStoreObj = getUpload(installationUrl, uploadId)
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
      updateUpload(installationUrl, uploadId, dataToUpdate)
    }).catch((error) => {
      if (requireEntriesPage && error.apiMessage === 'Page out of range requested.') {
        // Special case: can happen if entries have been deleted and the page we were on is no longer in range
        if (currentPagination && currentPagination.page !== 1) {
          // Rather than sending an update to all subscribers with an error, we first try
          // jumping to page 1 (will probably solve the problem)
          getUpload(installationUrl, uploadId).pagination.page = 1
          refreshUpload(installationUrl, uploadId)
          return
        }
      }
      updateUpload(installationUrl, uploadId, {error: error, isRefreshing: false, refreshOptions})
    })
  }

  /**
   * Use to nicely request a refresh of the upload store object.
   */
  function requestRefreshUpload(installationUrl, uploadId) {
    const uploadStoreObj = getUpload(installationUrl, uploadId)
    if (!uploadStoreObj.isRefreshing) {
      // Refresh is not already in progress
      refreshUpload(installationUrl, uploadId)
    }
  }

  async function initiateUploadRefreshIfNeeded(installationUrl, uploadId) {
    // Internal use: check if a refresh of the store is needed, and if so, initiate it.
    let uploadStoreObj = getUpload(installationUrl, uploadId)
    if (uploadStoreObj.isRefreshing) return // refresh already in progress
    if (uploadStoreObj.isProcessing) {
      // Upload is processing
      uploadStoreObj.isRefreshing = true // Signal start of a refresh
      await new Promise(resolve => setTimeout(resolve, 1000)) // wait one sec
      uploadStoreObj = getUpload(installationUrl, uploadId)
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
      refreshUpload(installationUrl, uploadId)
    } else {
      uploadStoreObj.isRefreshing = false
    }
  }

  /**
   * Gets an entry object from the store, creating it if it doesn't exist (in which case
   * an object with default values, mostly undefined or nulls, will be returned). Note, it
   * does not cause the store to fetch any data; for that, a subscription is necessary.
   */
  function getEntry(installationUrl, entryId) {
    if (!entryId) return undefined
    if (installationUrl !== apiBase) throw new Error('Fetching entries from external installations is not yet supported')
    let entryStoreObj = entryStore.current[entryId]
    if (!entryStoreObj) {
      // Creates an initial, empty entry store object.
      entryStoreObj = {
        installationUrl, // ReadOnly
        entryId: entryId, // ReadOnly
        uploadId: undefined, // ReadOnly - fetched by the store
        url: undefined, // ReadOnly - populated when uploadId fetched from the store
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
        isRefreshing: false,
        refreshOptions: null, // The options used at the last store refresh
        _subscriptions: [],

        // Convenience methods
        handleArchiveChanged: () => { handleArchiveChanged(installationUrl, entryId) },
        saveArchive: () => { return saveArchive(installationUrl, entryId) },
        reload: () => { requestRefreshEntry(installationUrl, entryId) }
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
  async function getEntryAsync(installationUrl, entryId, requireMetadata, requireArchive) {
    if (!entryId) return undefined
    const entryStoreObj = getEntry(installationUrl, entryId)
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
      subscribeToEntry(installationUrl, entryId, cb, requireMetadata, requireArchive)
    })
  }

  /**
   * Subscribes the callback cb to an entry, and returns a function to be called to unsubscribe.
   * Typically used in useEffect. The callback will be called when the store value changes.
   */
  function subscribeToEntry(installationUrl, entryId, cb, requireMetadata, requireArchive) {
    if (!entryId) return undefined
    if (requireMetadata === undefined || !(requireArchive === undefined || requireArchive === '*' || typeof requireArchive === 'object')) {
      throw Error('Store error: bad subscription parameter supplied')
    }
    const entryStoreObj = getEntry(installationUrl, entryId)
    addSubscription(entryStoreObj, cb, {requireMetadata, requireArchive})
    initiateEntryRefreshIfNeeded(installationUrl, entryId)
    return function unsubscriber() { removeSubscription(entryStore.current, entryId, cb) }
  }

  /**
   * Updates the store entry with the specified data and notifies all subscribers.
   */
  function updateEntry(installationUrl, entryId, dataToUpdate) {
    const oldStoreObj = getEntry(installationUrl, entryId)
    const newStoreObj = {...oldStoreObj, ...dataToUpdate}
    // Compute derived values not set by the refreshEntry method
    newStoreObj.url = newStoreObj.uploadId ? `${installationUrl}/uploads/${newStoreObj.uploadId}/archive/${entryId}` : undefined
    newStoreObj.exists = newStoreObj?.error?.name !== 'DoesNotExist'
    newStoreObj.archiveHasChanges = newStoreObj.archiveVersion !== newStoreObj.savedArchiveVersion

    // Update the store
    entryStore.current[entryId] = newStoreObj

    // Notify subscribers
    for (const subscription of [...newStoreObj._subscriptions]) {
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
    initiateEntryRefreshIfNeeded(installationUrl, entryId)
  }

  async function refreshEntry(installationUrl, entryId) {
    // Internal use: refresh an entry store obj with data from the API.
    let entryStoreObj = getEntry(installationUrl, entryId)
    let refreshOptions = entryOptions(entryStoreObj)
    let {requireMetadata, requireArchive} = refreshOptions
    if (!requireMetadata && !requireArchive) return
    // There's something that can be refreshed
    entryStoreObj.isRefreshing = true
    const dataToUpdate = {refreshOptions}
    try {
      if (requireMetadata) {
        const metadataApiData = await api.get(`/entries/${entryId}`, null, {returnRequest: true})
        const metadata = metadataApiData?.response?.data
        const uploadId = metadata?.upload_id
        const user = userRef.current
        const isWriter = user && metadata?.writers && metadata.writers.find(u => u.user_id === user.sub)
        const isEditableArchive = metadata && !metadata.published && metadata.quantities && metadata.quantities.includes('data')
        const editable = isWriter && isEditableArchive && selectedEntry.current === `${installationUrl}:${entryId}`
        const isProcessing = !!metadata?.process_running
        Object.assign(dataToUpdate, {metadataApiData, metadata, uploadId, editable, isProcessing, error: undefined})
        // Fetch the options again, in case some subscriptions were added while waiting for the api call
        entryStoreObj = getEntry(installationUrl, entryId)
        refreshOptions = entryOptions(entryStoreObj)
        requireArchive = refreshOptions.requireArchive
        dataToUpdate.refreshOptions.requireArchive = requireArchive
      }
      if (requireArchive) {
        const required = requireArchive === '*' ? '*' : {...requireArchive, 'resolve-inplace': false}
        const archiveApiData = await api.post(
          `/entries/${entryId}/archive/query`, {required}, {returnRequest: true, jsonResponse: true})
        const archive = archiveApiData?.response?.data?.archive
        if (archive) {
          archive.processing_logs = undefined
          if (archive.metadata?.upload_id) {
            dataToUpdate.uploadId = archive.metadata.upload_id
          }
        }
        dataToUpdate.archiveApiData = archiveApiData
        dataToUpdate.archive = archive
      }
    } catch (error) {
      dataToUpdate.error = error
    }
    dataToUpdate.isRefreshing = false
    updateEntry(installationUrl, entryId, dataToUpdate)
  }

  /**
   * Use to nicely request a refresh of the entry store object.
   */
  function requestRefreshEntry(installationUrl, entryId) {
    const entryStoreObj = getEntry(installationUrl, entryId)
    if (!entryStoreObj.isRefreshing) {
      // Refresh is not already in progress
      refreshEntry(installationUrl, entryId)
    }
  }

  async function initiateEntryRefreshIfNeeded(installationUrl, entryId) {
    // Internal use: check if a refresh of the store is needed, and if so, initiate it.
    let entryStoreObj = getEntry(installationUrl, entryId)
    if (entryStoreObj.isRefreshing) return // refresh already in progress
    if (entryStoreObj.isProcessing) {
      // Entry is processing
      entryStoreObj.isRefreshing = true // Signal start of a refresh
      await new Promise(resolve => setTimeout(resolve, 1000)) // wait one sec
      entryStoreObj = getEntry(installationUrl, entryId)
    }
    // Determine if a refresh is needed or not
    const {requireMetadata, requireArchive} = entryOptions(entryStoreObj)
    const lastRefreshSatisfiesOptions = entryRefreshSatisfiesOptions(entryStoreObj, requireMetadata, requireArchive)
    if (!entryStoreObj.error && (!lastRefreshSatisfiesOptions || entryStoreObj.isProcessing)) {
      // Need to fetch data from the api
      refreshEntry(installationUrl, entryId)
    } else {
      entryStoreObj.isRefreshing = false
    }
  }

  /**
   * Used to save the archive and trigger a store refresh.
   */
  function saveArchive(installationUrl, entryId) {
    const {uploadId, metadata, archive, archiveVersion} = getEntry(installationUrl, entryId)
    const {mainfile} = metadata
    if (uploadId) {
      const separatorIndex = mainfile.lastIndexOf('/')
      const path = mainfile.substring(0, separatorIndex + 1)
      const fileName = mainfile.substring(separatorIndex + 1)
      const newArchive = {...archive}
      delete newArchive.metadata
      delete newArchive.results
      delete newArchive.processing_logs

      const config = {}
      let stringifiedArchive
      if (fileName.endsWith('yaml') || fileName.endsWith('yml')) {
        config.headers = {
          'Content-Type': 'application/yaml'
        }
        stringifiedArchive = YAML.stringify(newArchive)
      }
      return new Promise((resolve, reject) => {
        api.put(`/uploads/${uploadId}/raw/${path}?file_name=${fileName}&wait_for_processing=true&entry_hash=${archive.metadata.entry_hash}`, stringifiedArchive || newArchive, config)
          .then(response => {
            requestRefreshEntry(installationUrl, entryId)
          })
          .catch(error => {
            if (error?.status === 409) {
              reject(error)
            } else {
              raiseError(error)
            }
          })
        updateEntry(installationUrl, entryId, {savedArchiveVersion: archiveVersion})
      })
    }
  }

  /**
   * Call to signal that the archive has been manually edited.
   */
  function handleArchiveChanged(installationUrl, entryId) {
    const {archiveVersion} = getEntry(installationUrl, entryId)
    updateEntry(installationUrl, entryId, {archiveVersion: archiveVersion + 1})
  }

  /**
   * Gets a parsed metainfo object from the store, asynchronously.
   * The url can be a string or a parsed Url object. If it is empty, null will be returned.
   * Note, this method always returns the whole data object. I.e. if the url specifies a
   * particular section definition, we return the metainfo data object which *contains* this
   * definition, not just the definition itself.
   */
  async function getMetainfoAsync(url) {
    if (!url) return null
    // Get metainfoData
    const metainfoBaseUrl = getMetainfoBaseUrl(url)
    let metainfoData = metainfoDataStore.current[metainfoBaseUrl]
    if (!metainfoData) {
      // Not found in cache. Fetch data object.
      const parsedMetainfoBaseurl = parseNomadUrl(metainfoBaseUrl)
      metainfoData = await fetchMetainfoData(parsedMetainfoBaseurl)
      metainfoData._url = metainfoBaseUrl
      metainfoData._parsedUrl = parsedMetainfoBaseurl
      metainfoData._subscriptions = []
      metainfoDataStore.current[metainfoBaseUrl] = metainfoData
    }
    // If needed, create metainfo (which will also initiate parsing)
    if (!metainfoData._metainfo) {
      const parent = metainfoBaseUrl === systemMetainfoUrl ? undefined : await getMetainfoAsync(systemMetainfoUrl)
      metainfoData._metainfo = new Metainfo(parent, metainfoData, getMetainfoAsync, externalInheritanceCache.current)
    }
    // Returned object after parsing is completed.
    return await metainfoData._metainfo._result
  }

  /**
   * Gets the parsed metainfo definition specified by the url, asynchronously.
   * The url can be a string or a parsed Url object. If it is empty, null will be returned.
   */
  async function getMetainfoDefAsync(url) {
    if (!url) return null
    const parsedUrl = parseNomadUrl(url)
    const metainfo = await getMetainfoAsync(parsedUrl)
    if (parsedUrl.qualifiedName) {
      return metainfo.getDefByQualifiedName(parsedUrl.qualifiedName)
    } else if (parsedUrl.path) {
      return metainfo.getDefByPath(parsedUrl.path)
    }
    throw new Error(`Expected metainfo definition url, got: ${url}`)
  }

  /**
   * Subscribes the callback cb to a metainfo, and returns a function to be called to unsubscribe.
   * The main purpose of the subscriptions is to tell the store that it needs to keep this
   * metainfo object in the cache. If a metainfo object is not subscribed to by anyone, it will
   * eventually be removed from the cache (with the exception of the system metainfo, which is
   * always kept). If the url is empty, or refers to the system metainfo, we do nothing and
   * return null.
   */
  function subscribeToMetainfo(url, cb) {
    if (!url) return null
    const metainfoBaseUrl = getMetainfoBaseUrl(url)
    if (metainfoBaseUrl === systemMetainfoUrl) return null
    const metainfoData = metainfoDataStore.current[metainfoBaseUrl]
    for (const url of metainfoData._metainfo._dependencies) {
      const dependency = metainfoDataStore.current[url]
      if (dependency) {
        addSubscription(dependency, cb)
      }
    }
    return function unsubscriber() {
      const metainfoData = metainfoDataStore.current[metainfoBaseUrl]
      for (const url of metainfoData._metainfo._dependencies) {
        removeSubscription(metainfoDataStore.current, url, cb)
      }
    }
  }

  function getMetainfoBaseUrl(url) {
    // Internal use: computes the *metainfo base url* (string), given an absolute metainfo url.
    // The base url is an absolute url, identifying the metainfo data obj to parse
    // (to get the Metainfo object), and is also used as the cache key to store the resulting
    // Metainfo object.
    if ((url.url || url) === systemMetainfoUrl) return systemMetainfoUrl
    url = parseNomadUrl(url)
    if (url.type !== refType.metainfo && url.type !== refType.archive) throw new Error(`Cannot get metainfo: bad url type (${url})`)
    if (!url.isResolved) throw new Error(`Url not resolved: ${url}`)

    if (url.versionHash) {
      // TODO
      throw new Error('Fetching frozen schemas not yet implemented')
    } else if (url.entryId) {
      return createEntryUrl(url.installationUrl, url.uploadId, url.entryId, '/definitions')
    } else if (url.qualifiedName) {
      return systemMetainfoUrl
    }
  }

  async function fetchMetainfoData(parsedMetainfoBaseUrl) {
    // Internal use: fetches the metainfo data
    if (parsedMetainfoBaseUrl.url === systemMetainfoUrl) return currentSystemMetainfoData
    if (parsedMetainfoBaseUrl.versionHash) {
      // TODO
      throw new Error('Fetching frozen schemas not yet implemented')
    } else if (parsedMetainfoBaseUrl.entryId) {
      const entryStoreObj = await getEntryAsync(
        parsedMetainfoBaseUrl.installationUrl, parsedMetainfoBaseUrl.entryId, false, metainfoArchiveFilter)
      if (entryStoreObj.error) {
        throw new Error(`Error fetching entry ${parsedMetainfoBaseUrl.entryId}: ${entryStoreObj.error}`)
      } else if (!entryStoreObj.archive?.definitions) {
        throw new Error(`Entry ${parsedMetainfoBaseUrl.entryId} does not contain metainfo definitions`)
      }
      return {
        definitions: JSON.parse(JSON.stringify(entryStoreObj.archive.definitions)),
        metadata: JSON.parse(JSON.stringify(entryStoreObj.archive.metadata))
      }
    }
  }

  /**
   * Gets all the inheriting sections of the provided section definition. That is: all inheriting
   * sections *currently known to the store*. The result is returned as a list of definitions.
   */
  function getAllInheritingSections(definition) {
    const rv = []
    // Add subclasses (i.e. *directly* inheriting from this section)
    const url = getUrlFromDefinition(definition)
    rv.push(...(definition._allInternalInheritingSections || []))
    rv.push(...(externalInheritanceCache.current[url] || []))
    // Add recursively (sub-sub classes, sub-sub-sub classes etc.)
    const recursive = []
    for (const inheritingSection of rv) {
      recursive.push(...getAllInheritingSections(inheritingSection))
    }
    rv.push(...recursive)
    return rv
  }

  const contextValue = {
    getUpload,
    getUploadAsync,
    subscribeToUpload,
    updateUpload,
    requestRefreshUpload,
    getEntry,
    getEntryAsync,
    subscribeToEntry,
    selectedEntry,
    getMetainfoAsync,
    getMetainfoDefAsync,
    subscribeToMetainfo,
    getAllInheritingSections
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

/**
 * React function for getting an entry from the store with certain subscription options.
 * If the store already has fetched the entry and all the required data, it will be returned
 * immediately from the cache. If the store has loaded the entry, but it has not yet fetched
 * all data required by the subscription options, the function will return a copy of the store
 * object without the following fields: metadata, metadataApiData, archive, archiveApiData.
 * This is to signal to the component that all data requested is not yet available in the store.
 *
 * Subscription options:
 * @param {*} requireMetadata If we require the store to fetch the entry metadata
 * @param {*} requireArchive If we require the store to fetch the archive, and if so, what
 *    parts of it. Should be one of:
 *      a) undefined (no archive data required)
 *      b) '*' (load entire archive), or
 *      c) an object specifying a simple archive data filter (see the doc for mergeArchiveFilter).
 */
export function useEntryStoreObj(installationUrl, entryId, requireMetadata, requireArchive) {
  const dataStore = useDataStore()
  const [entryStoreObj, setEntryStoreObj] = useState(
    () => installationUrl && entryId
      ? filteredEntryStoreObj(dataStore.getEntry(installationUrl, entryId), requireMetadata, requireArchive)
      : null)

  const onEntryStoreUpdated = useCallback((oldStoreObj, newStoreObj) => {
    setEntryStoreObj(filteredEntryStoreObj(newStoreObj, requireMetadata, requireArchive))
  }, [setEntryStoreObj, requireMetadata, requireArchive])

  useEffect(() => {
    if (installationUrl && entryId) {
      return dataStore.subscribeToEntry(installationUrl, entryId, onEntryStoreUpdated, requireMetadata, requireArchive)
    }
  }, [installationUrl, entryId, requireMetadata, requireArchive, dataStore, onEntryStoreUpdated])

  return entryStoreObj
}

/**
 * Misc internal helper funcions
 */

function filteredEntryStoreObj(entryStoreObj, requireMetadata, requireArchive) {
  // Returns a filtered entry store obj if all data is not yet available.
  if (!entryRefreshSatisfiesOptions(entryStoreObj, requireMetadata, requireArchive)) {
    const rv = {...entryStoreObj}
    delete rv.metadata
    delete rv.metadataApiData
    delete rv.archive
    delete rv.archiveApiData
    return rv
  }
  return entryStoreObj
}

function uploadOptions(uploadStoreObj) {
  // Internal use: Determine the options of the upload, based on the subscriptions
  let requireUpload = false
  let requireEntriesPage = false
  for (const subscription of uploadStoreObj._subscriptions) {
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

function entryOptions(entryStoreObj) {
  // Internal use: Determine the options of the entry, based on the subscriptions
  let requireMetadata = false
  let requireArchive
  for (const subscription of entryStoreObj._subscriptions) {
    requireMetadata = requireMetadata || subscription.requireMetadata
    requireArchive = mergeArchiveFilter(requireArchive, subscription.requireArchive)
  }
  return {requireMetadata, requireArchive}
}

function entryRefreshSatisfiesOptions(entryStoreObj, requireMetadata, requireArchive) {
  // Internal use: Determine if the last store refresh satisfies the given options.
  const refreshOptions = entryStoreObj?.refreshOptions
  if (refreshOptions) {
    return (refreshOptions.requireMetadata || !requireMetadata) && archiveFilterExtends(refreshOptions.requireArchive, requireArchive)
  }
  return false
}

/**
 * Merges two archive data filters, filter1 and filter2, returning a filter which will fetch all data
 * fetched by either filter1 or filter2.
 * Note, this method only handles certain simple types of filters, used by the gui.
 * Each filter can either have the value
 *  1) undefined (meaning include nothing)
 *  2) '*' (meaning include everything, but do not add resolved data)
 *  3) 'include-resolved' (include everything and add resolved data)
 *  4) an object with subkeys (will be treated recursively)
 * The `resolve-inplace` option should not be specified (it will be implicitly set to false).
 */
function mergeArchiveFilter(filter1, filter2, level = 0) {
  if (filter1 === filter2) return filter1
  if (filter1 === undefined) return filter2
  if (filter2 === undefined) return filter1
  if (filter1 === 'include-resolved' || filter2 === 'include-resolved') return 'include-resolved'
  const obj1 = typeof filter1 === 'object'
  const obj2 = typeof filter2 === 'object'
  if (filter1 === '*' && obj2) return level !== 0 && hasIncludeResolved(filter2) ? 'include-resolved' : '*'
  if (filter2 === '*' && obj1) return level !== 0 && hasIncludeResolved(filter1) ? 'include-resolved' : '*'
  if (!obj1 || !obj2) throw new Error('Cannot merge filters: bad values')

  // Two objects. Inspect recursively.
  const rv = {}
  for (const key in {...filter1, ...filter2}) {
    const value1 = filter1[key]
    const value2 = filter2[key]
    rv[key] = mergeArchiveFilter(value1, value2, level + 1)
  }
  return rv
}

/**
 * True if the archive filter filter1 *extends* the filter filter2, i.e. if all data passing
 * filter2 will also pass filter1.
 * Note, this method only handles certain simple types of filters, used by the gui.
 */
function archiveFilterExtends(filter1, filter2, level = 0) {
  if (filter1 === filter2) return true
  if (filter2 === undefined) return true
  if (filter1 === undefined) return false
  if (filter1 === 'include-resolved') return true
  if (filter2 === 'include-resolved') return false
  const obj1 = typeof filter1 === 'object'
  const obj2 = typeof filter2 === 'object'
  if (filter1 === '*' && obj2) return level === 0 || !hasIncludeResolved(filter2)
  if (filter2 === '*' && obj1) return false
  if (!obj1 || !obj2) throw new Error('Cannot compare filters: bad values')
  // Two objects. Inspect recursively.
  for (const key in filter2) {
    const value1 = filter1[key]
    const value2 = filter2[key]
    if (!archiveFilterExtends(value1, value2, level + 1)) return false
  }
  return true
}

function hasIncludeResolved(f) {
  for (const key in f) {
    const v = f[key]
    if (v === 'include-resolved') return true
    if (typeof v === 'object') {
      if (hasIncludeResolved(v)) return true
    }
  }
  return false
}
