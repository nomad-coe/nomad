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
import { useApi } from '../api'
import { useErrors } from '../errors'

const entryContext = React.createContext()

export const useEntryContext = () => {
  return useContext(entryContext)
}

const EntryContext = React.memo(function EntryContext({entryId, children}) {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [metadataApiData, setMetadataApiData] = useState()
  const metadata = useMemo(() => metadataApiData?.response?.data, [metadataApiData])
  const [archiveRequired, setArchiveRequired] = useState(false)
  const [archiveApiData, setArchiveApiData] = useState()
  const [archiveVersion, setArchiveVersion] = useState(0)
  const [savedArchiveVersion, setSavedArchiveVersion] = useState(0)

  const archive = useMemo(() => {
    const archive = archiveApiData?.response?.data?.archive
    if (archive) {
      archive.processing_logs = undefined
    }
    return archive
  }, [archiveApiData])
  const [exists, setExists] = useState(true)

  const fetchMetadata = useCallback(() => {
    api.get(`/entries/${entryId}`, null, {returnRequest: true}).then(setMetadataApiData).catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, entryId, setMetadataApiData, setExists])

  useEffect(() => {
    fetchMetadata()
  }, [fetchMetadata])

  useEffect(() => {
    if (archiveRequired && metadata && exists) {
      api.post(
        `/entries/${entryId}/archive/query`,
        {
          required: archiveRequired
        }, {returnRequest: true, jsonResponse: true}
      )
        .then(setArchiveApiData)
        .catch(error => {
          if (error.name === 'DoesNotExist') {
            setExists(false)
          } else {
            raiseError(error)
          }
        })
    }
  }, [archiveRequired, metadata, exists, setExists, setArchiveApiData, api, entryId, raiseError])

  const requireArchive = useCallback((required) => {
    setArchiveRequired(required || '*')
  }, [setArchiveRequired])

  const reload = useCallback(() => {
    fetchMetadata()
  }, [fetchMetadata])

  const saveArchive = useCallback((callback) => {
    const uploadId = metadata.upload_id
    const {mainfile} = metadata
    if (uploadId) {
      const separatorIndex = mainfile.lastIndexOf('/')
      const path = mainfile.substring(0, separatorIndex + 1)
      const fileName = mainfile.substring(separatorIndex + 1)
      const newArchive = {...archive}
      delete newArchive.metadata
      delete newArchive.results
      delete newArchive.processing_logs
      delete newArchive.resources
      api.put(`/uploads/${uploadId}/raw/${path}?file_name=${fileName}&wait_for_processing=true`, newArchive)
        .then(response => {
          // TODO handle processing errors
          reload()
        })
        .catch(raiseError)
    }
    setSavedArchiveVersion(archiveVersion)
  }, [
    setSavedArchiveVersion, archiveVersion, api, reload,
    raiseError, archive, metadata
  ])

  const handleArchiveChanged = useCallback(() => {
    setArchiveVersion(value => value + 1)
  }, [setArchiveVersion])

  const value = useMemo(() => {
    const editable = metadata && !metadata.published && (metadata.quantities.find(quantity => {
      return quantity === 'data'
    }) !== null)
    return {
      entryId: entryId,
      uploadId: metadata?.upload_id,
      exists: exists,
      metadata: metadata,
      metadataApiData: metadataApiData,
      archive: archive,
      archiveApiData: archiveApiData,
      editable: editable,
      requireArchive: requireArchive,
      reload: reload,
      saveArchive: saveArchive,
      handleArchiveChanged: handleArchiveChanged,
      archiveHasChanges: archiveVersion !== savedArchiveVersion,
      archiveVersion: archiveVersion
    }
  }, [
    entryId, exists, metadata, archive, metadataApiData, archiveApiData,
    savedArchiveVersion, archiveVersion, requireArchive, reload, handleArchiveChanged,
    saveArchive
  ])
  return <entryContext.Provider value={value}>
    {children}
  </entryContext.Provider>
})
EntryContext.propTypes = {
  entryId: PropTypes.string.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default EntryContext
