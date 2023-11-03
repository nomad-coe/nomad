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
import React, { useState, useMemo, useEffect } from 'react'
import PropTypes from 'prop-types'
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
import { MaterialRowActions } from '../material/MaterialDetails'
import { pluralize, formatInteger } from '../../utils'
import { isEmpty } from 'lodash'
import { useSearchContext } from './SearchContext'

/**
 * Displays the list of search results.
 */
const SearchResults = React.memo((props) => {
  const {noAction, onSelectedChanged, defaultUncollapsedEntryID, title, 'data-testid': testID, PaperProps, ...otherProps} = props
  const {columns, resource, rows, useResults, useApiQuery} = useSearchContext()
  const {data, pagination, setPagination} = useResults()
  const apiQuery = useApiQuery()
  const [selected, setSelected] = useState(new Set())

  useEffect(() => {
    if (onSelectedChanged) {
      onSelectedChanged(selected)
    }
  }, [onSelectedChanged, selected])

  const query = useMemo(() => {
    if (selected === 'all') {
      return apiQuery
    }
    return {entry_id: [...selected]}
  }, [selected, apiQuery])

  if (isEmpty(columns)) {
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

  // Select components based on the targeted resource
  let details
  let actions
  let buttons
  if (resource === "entries") {
    details = rows?.details?.render || EntryDetails
    actions = rows?.actions?.render || EntryRowActions
    if (!noAction) buttons = <EntryDownloadButton tooltip="Download files" query={query} />
  } else if (resource === "materials") {
    actions = MaterialRowActions
  }

  return <Paper data-testid={testID} {...PaperProps}>
    <Datatable
      data={data}
      pagination={pagination}
      onPaginationChanged={setPagination}
      columns={columns?.options && Object.values(columns.options)}
      shownColumns={columns?.selected}
      selected={rows?.selection?.enabled ? selected : undefined}
      getId={option => option.entry_id}
      onSelectedChanged={rows?.selection?.enabled ? setSelected : undefined}
      {...otherProps}
    >
      <DatatableToolbar title={`${formatInteger(data.length)}/${pluralize(title || 'result', pagination.total, true, true, title ? undefined : 'search')}`}>
        {rows?.selection?.enabled &&
          <DatatableToolbarActions selection>
            {buttons}
          </DatatableToolbarActions>
        }
      </DatatableToolbar>
      <DatatableTable
        actions={rows?.actions?.enabled ? actions : undefined}
        details={rows?.details?.enabled ? details : undefined}
        defaultUncollapsedRow={defaultUncollapsedEntryID && data.find(row => row.entry_id === defaultUncollapsedEntryID)}
      >
        <DatatableLoadMorePagination color="primary">load more</DatatableLoadMorePagination>
      </DatatableTable>
    </Datatable>
  </Paper>
})

SearchResults.propTypes = {
  noAction: PropTypes.bool,
  PaperProps: PropTypes.object,
  onSelectedChanged: PropTypes.func,
  defaultUncollapsedEntryID: PropTypes.string,
  title: PropTypes.string,
  'data-testid': PropTypes.string
}

SearchResults.defaultProps = {
  'data-testid': 'search-results'
}

export default SearchResults
