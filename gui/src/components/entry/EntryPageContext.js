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

import React, { useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import { useDataStore, useEntryStoreObj } from '../DataStore'
import { apiBase } from '../../config'

const entryPageContext = React.createContext()

export const useEntryPageContext = (requireArchive) => {
  const entryId = useContext(entryPageContext)
  return useEntryStoreObj(apiBase, entryId, true, requireArchive)
}

const EntryPageContext = React.memo(function EntryContext({entryId, children}) {
  const {selectedEntry} = useDataStore()

  useEffect(() => {
    // Inform the Store of the selected entry
    selectedEntry.current = `${apiBase}:${entryId}`
    return () => { selectedEntry.current = undefined }
  }, [selectedEntry, entryId])

  return <entryPageContext.Provider value={entryId}>
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
