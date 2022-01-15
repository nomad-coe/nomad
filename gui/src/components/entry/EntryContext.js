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
  const [loadArchive, setLoadArchive] = useState(false)
  const [archiveApiData, setArchiveApiData] = useState()

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
    if (loadArchive && metadata && exists) {
      api.get(`/entries/${entryId}/archive`, null, {returnRequest: true})
        .then(setArchiveApiData)
        .catch(error => {
          if (error.name === 'DoesNotExist') {
            setExists(false)
          } else {
            raiseError(error)
          }
        })
    }
  }, [loadArchive, metadata, exists, setExists, setArchiveApiData, api, entryId, raiseError])

  const requireArchive = useCallback(() => {
    setLoadArchive(true)
  }, [setLoadArchive])

  const reload = useCallback(() => {
    fetchMetadata()
  }, [fetchMetadata])

  const value = useMemo(() => ({
    entryId: entryId,
    exists: exists,
    metadata: metadata,
    metadataApiData: metadataApiData,
    archive: archive,
    archiveApiData: archiveApiData,
    editable: metadata?.published,
    requireArchive: requireArchive,
    reload: reload
  }), [entryId, exists, metadata, archive, metadataApiData, archiveApiData, requireArchive, reload])
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
