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

import React, { useContext, useState, useEffect, useRef, useMemo } from 'react'
import PropTypes from 'prop-types'
import { cloneDeep } from 'lodash'
import { useDataStore, useEntryStoreObj } from '../DataStore'
import { apiBase, ui } from '../../config'

const entryPageContext = React.createContext()

/**
 * Hook for fetching data from the current EntryPageContext.
 *
 * @param {*} requireArchive Optional query filter
 * @param {*} update Whether to keep updating the data for the original entry
 *   if changes beyond the first successfull archive load are made in the store.
 *   Sometimes the first version of the data should be kept to avoid unnecessary
 *   re-renders and layout inconsistencies. Note that if the entry id changes,
 *   the update is forced.
 * @returns
 */
export const useEntryPageContext = (requireArchive, update = true) => {
  const {entryId, overview} = useContext(entryPageContext) || {}
  const entryData = useEntryStoreObj(apiBase, entryId, true, requireArchive)
  const [data, setData] = useState(entryData)
  const oldEntryId = useRef()
  const completed = useRef(false)

  // Get the overview config
  const finalOverview = useMemo(() => {
    const tmpOverview = overview || ui?.entry_context?.overview
    if (!tmpOverview) return {}
    const finalOverview = cloneDeep(tmpOverview)
    const options = finalOverview.include
      .filter(key => !finalOverview?.exclude?.includes(key))
      .map(key => ({key, ...finalOverview.options[key]}))
    return {options}
  }, [overview])

  // This effect controls how the data returned by this hook is synchronized
  // with the data coming from the store. If update = true, the data is always
  // synchronized. If update = false, it is only synchronized until the archive
  // data is fully set, or the entry changes.
  useEffect(() => {
    const newEntryId = entryData?.entryId
    const isNew = newEntryId !== oldEntryId.current
    if (update) {
      setData(entryData)
    } else if (isNew || !completed.current) {
      setData(entryData)
      if (isNew) {
        completed.current = false
      }
      if (entryData?.archive) {
        completed.current = true
      }
    }
    oldEntryId.current = newEntryId
  }, [entryData, update])

  // The final data is memoized in order to avoid unwanted rerenders
  const context = useMemo(() => {
    return {...data, overview: finalOverview}
  }, [data, finalOverview])

  return context
}

const EntryPageContext = React.memo(({entryId, overview, children}) => {
  const dataStore = useDataStore()
  dataStore.resetIfNeeded(entryId)

  useEffect(() => {
    // Inform the Store of the selected entry
    dataStore.selectedEntry.current = `${apiBase}:${entryId}`
    return () => { dataStore.selectedEntry.current = undefined }
  }, [dataStore, entryId])

  const value = useMemo(() => {
    return {
      entryId,
      overview
    }
  }, [entryId, overview])

  return <entryPageContext.Provider value={value}>
    {children}
  </entryPageContext.Provider>
})
EntryPageContext.propTypes = {
  entryId: PropTypes.string.isRequired,
  overview: PropTypes.object,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default EntryPageContext
