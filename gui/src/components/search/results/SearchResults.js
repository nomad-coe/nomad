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
import React from 'react'
import {
  Paper, Typography
} from '@material-ui/core'
import SearchResultsMaterials from './SearchResultsMaterials'
import SearchResultsEntries from './SearchResultsEntries'
import { useScrollResults, useSearchContext } from '../SearchContext'

const orderByMap = {
  'entries': 'upload_create_time',
  'materials': 'chemical_formula_hill'
}

/**
 * Displays the list of search results
 */
const SearchResults = React.memo(function SearchResults() {
  const {resource} = useSearchContext()
  const {data, pagination, setPagination} = useScrollResults({
    page_size: 20, order_by: orderByMap[resource]
  })

  if (pagination.total === 0) {
    return <Typography>no results</Typography>
  }

  if (!pagination.total) {
    return <Typography>searching ...</Typography>
  }

  const Component = resource === 'materials' ? SearchResultsMaterials : SearchResultsEntries
  return <Paper>
    <Component
      data={data}
      pagination={pagination} onPaginationChanged={setPagination}
    />
  </Paper>
})

export default SearchResults
