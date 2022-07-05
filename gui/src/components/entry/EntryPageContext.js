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
import { useDataStore } from '../DataStore'

const entryPageContext = React.createContext()

export const useEntryPageContext = () => {
  return useContext(entryPageContext)
}

const EntryPageContext = React.memo(function EntryContext({entryId, children}) {
  const dataStore = useDataStore()
  const [entryStoreObj, setEntryStoreObj] = useState(dataStore.getEntry(entryId))

  const onEntryStoreUpdated = useCallback((oldStoreObj, newStoreObj) => {
    setEntryStoreObj(newStoreObj)
  }, [setEntryStoreObj])

  useEffect(() => {
    return dataStore.subscribeToEntry(entryId, onEntryStoreUpdated, true, true, true)
  }, [dataStore, entryId, onEntryStoreUpdated])

  const contextValue = useMemo(() => { return entryStoreObj }, [entryStoreObj])

  return <entryPageContext.Provider value={contextValue}>
    {children}
  </entryPageContext.Provider>
})
EntryPageContext.propTypes = {
  entryId: PropTypes.string.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default EntryPageContext
