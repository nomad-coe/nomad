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

import React, { useState, useMemo, useEffect, useContext } from 'react'
import PropTypes from 'prop-types'
import { useApi } from '../api'
import { useDataStore, useEntryStoreObj } from '../DataStore'
import { ui, apiBase} from '../../config'

/**
 * Hook for retrieving read-only archive information for the entry present in
 * the active EntryContext.
 */
export const useArchive = (required) => {
  const { entryId } = useEntryContext()
  const { api } = useApi()
  const [response, setResponse] = useState({})

  // When entryId changes, resets the currently held response and updates it
  // once new data arrives
  useEffect(() => {
    setResponse({})
    if (!entryId) return
    api
      .post(`entries/${entryId}/archive/query`, {required}, {returnRequest: true})
      .then(response => setResponse({data: response?.response?.data?.archive, response: response}))
  }, [api, required, entryId])

  return response
}

/**
 * Hook for retrieving read-only index information for the entry present in the
 * active EntryContext.
 */
export const useIndex = () => {
  const { entryId } = useEntryContext()
  const { api } = useApi()
  const [response, setResponse] = useState({})

  // When entryId changes, resets the currently held response and updates it
  // once new data arrives
  useEffect(() => {
    setResponse({})
    api
      .get(`entries/${entryId}`, null, {returnRequest: true})
      .then(response => setResponse({data: response?.response?.data, response: response}))
  }, [api, entryId])

  return response
}

/**
 * Hook for fetching the "full" entry data. If the entry is editable, the user
 * can also mutate the archive through this hook.
 *
 * @param {*} requireArchive Optional query filter
 */
export const useEntryStore = (requireArchive) => {
  const {entryId} = useContext(entryContext) || {}
  const entryStoreObj = useEntryStoreObj(apiBase, entryId, true, requireArchive)
  return entryStoreObj
}

/**
 * Hook for fetching the current entry context.
 */
export const useEntryContext = () => {
  return useContext(entryContext)
}

/**
 * Context that manages information about a specific entry that will be targeted
 * by different components.
 */
const entryContext = React.createContext()
export const EntryContext = React.memo(({entryId, cards, children}) => {
  const dataStore = useDataStore()
  dataStore.resetIfNeeded(entryId)

  // Inform the Store of the selected entry
  useEffect(() => {
    dataStore.selectedEntry.current = `${apiBase}:${entryId}`
    return () => { dataStore.selectedEntry.current = undefined }
  }, [dataStore, entryId])

  // Get the overview config
  const finalCards = useMemo(() => {
    const finalCards = cards || ui?.entry?.cards
    return finalCards || undefined
  }, [cards])

  const value = useMemo(() => ({entryId, cards: finalCards}), [entryId, finalCards])
  return <entryContext.Provider value={value}>
    {children}
  </entryContext.Provider>
})
EntryContext.propTypes = {
  entryId: PropTypes.string,
  cards: PropTypes.object,
  children: PropTypes.node
}
