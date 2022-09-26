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
import React, { useState, useMemo } from 'react'
import { Paper, Typography } from '@material-ui/core'
import { Alert } from '@material-ui/lab'
import {
  Datatable,
  DatatableLoadMorePagination,
  DatatableTable,
  DatatableToolbar,
  DatatableToolbarActions
} from '../datatable/Datatable'
import EntryDownloadButton from '../entry/EntryDownloadButton'
import EntryDetails, { EntryRowActions } from '../entry/EntryDetails'
import { pluralize, formatInteger } from '../../utils'
import { useSearchContext } from './SearchContext'

/**
 * Displays the list of search results.
 */
const SearchResults = React.memo(function SearchResults() {
  const {columns, useResults, useQuery} = useSearchContext()
  const {data, pagination, setPagination} = useResults()
  const searchQuery = useQuery()
  const [selected, setSelected] = useState([])

  const query = useMemo(() => {
    if (selected === 'all') {
      return searchQuery
    }
    return {entry_id: selected.map(data => data.entry_id)}
  }, [selected, searchQuery])

  if (!columns) {
    return <Alert severity="warning">
      No search columns defined within this search context. Ensure that all GUI artifacts are created.
    </Alert>
  }

  if (pagination.total === 0) {
    return <Typography>no results</Typography>
  }

  if (!pagination.total) {
    return <Typography>searching ...</Typography>
  }

  return <Paper>
    <Datatable
      data={data}
      pagination={pagination}
      onPaginationChanged={setPagination}
      columns={columns?.options}
      shownColumns={columns?.enable}
      selected={selected}
      onSelectedChanged={setSelected}
    >
      <DatatableToolbar title={`${formatInteger(data.length)}/${pluralize('result', pagination.total, true, true, 'search')}`}>
        <DatatableToolbarActions selection>
          <EntryDownloadButton tooltip="Download files" query={query} />
        </DatatableToolbarActions>
      </DatatableToolbar>
      <DatatableTable actions={EntryRowActions} details={EntryDetails}>
        <DatatableLoadMorePagination color="primary">load more</DatatableLoadMorePagination>
      </DatatableTable>
    </Datatable>
  </Paper>
})

export default SearchResults
