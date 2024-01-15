/* eslint-disable quotes */
/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
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

import React, { memo, useMemo } from 'react'
import PropTypes from 'prop-types'
import { PropertyCard } from './PropertyCard'
import { SearchContext, useSearchContext } from "../../search/SearchContext"
import SearchResults from "../../search/SearchResults"
import { cloneDeep } from "lodash"
import { ui } from "../../../config"

const context = cloneDeep(ui?.apps?.options?.entries)
const shownColumns = ['entry_name', 'entry_type', 'entry_create_time']

const options = {
  'entry_name': {'label': 'Name', 'align': 'left'},
  'entry_type': {'label': 'Entry type', 'align': 'left'},
  'entry_create_time': {'label': 'Entry creation time', 'align': 'left'},
  'upload_name': {'label': 'Upload name', 'align': 'left'},
  'upload_id': {'label': 'Upload id', 'align': 'left'},
  'upload_create_time': {'label': 'Upload time', 'align': 'left'},
  'authors': {'label': 'Authors', 'align': 'left', 'sortable': false},
  'results.eln.lab_ids': {'label': 'IDs'}
}

Object.keys(options).forEach(key => { options[key].key = key })
context.columns.options = options
context.columns.selected = shownColumns
context.rows.details = { enabled: false }
context.rows.actions = { enabled: true }

const HistoryCard = memo(() => {
  const { useResults } = useSearchContext()
  const { data } = useResults()

  if (!data || data.length === 0) {
    return null
  }

  return (
    <PropertyCard title='History'>
      <SearchResults
        title='activity'
        multiSelect={false}
        PaperProps={{elevation: 0}}
      />
    </PropertyCard>
  )
})

const SampleHistoryUsingCard = memo(({ index }) => {
  const filtersLocked = useMemo(() => ({
    'entry_references.target_entry_id': [index.entry_id],
    'section_defs.definition_qualified_name': ['nomad.datamodel.metainfo.basesections.Activity']
  }), [index.entry_id])

  if (!index.section_defs?.some(def => def?.definition_qualified_name === 'nomad.datamodel.metainfo.basesections.Entity')) {
    return null
  }

  return (
    <SearchContext
      resource={context?.resource}
      initialPagination={context?.pagination}
      initialColumns={context?.columns}
      initialRows={context?.rows}
      initialFilterMenus={context?.filter_menus}
      initialFiltersLocked={filtersLocked}
      initialSearchSyntaxes={context?.search_syntaxes}
    >
      <HistoryCard />
    </SearchContext>
  )
})

SampleHistoryUsingCard.propTypes = {
  index: PropTypes.shape({
    entry_id: PropTypes.string.isRequired,
    section_defs: PropTypes.arrayOf(PropTypes.shape({
      definition_qualified_name: PropTypes.string.isRequired
    })).isRequired,
    entry_references: PropTypes.arrayOf(PropTypes.shape({
      target_entry_id: PropTypes.string.isRequired,
      target_path: PropTypes.string.isRequired
    }))
  }).isRequired
}

export default SampleHistoryUsingCard
