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
import React, { useState, useMemo, useEffect, useCallback } from 'react'
import PropTypes from 'prop-types'
import jmespath from 'jmespath'
import { Paper, Typography, Tooltip, IconButton } from '@material-ui/core'
import { Alert } from '@material-ui/lab'
import LaunchIcon from '@material-ui/icons/Launch'
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
import { pluralize, formatInteger, parseJMESPath, filterOptions } from '../../utils'
import { isEmpty, isArray } from 'lodash'
import { useSearchContext } from './SearchContext'

/**
 * Used to retrieve an URL link from the row metadata and display a link icon to
 * that resource.
 */
export const ActionURL = React.memo(({action, data}) => {
  const {path} = parseJMESPath(action.path)
  let href = jmespath.search(data, path)
  href = isArray(href) ? href[0] : href
  return <Tooltip title={action.description || ''}>
    <IconButton href={href} target="_blank"><LaunchIcon/></IconButton>
  </Tooltip>
})
ActionURL.propTypes = {
  action: PropTypes.object.isRequired, // Action configuration from app config
  data: PropTypes.object.isRequired // ES index data
}

/**
 * Displays the list of search results.
 */
export const SearchResults = React.memo(({
  columns,
  resource,
  rows,
  useResults,
  useApiQuery,
  noAction,
  onSelectedChanged,
  defaultUncollapsedEntryID,
  title,
  'data-testid': testID,
  PaperProps,
  ...otherProps
}) => {
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

  const actions = useCallback((data) => {
    const actionComponents = []
    for (const [key, value] of Object.entries(filterOptions(rows.actions))) {
      const component = {
        'url': <ActionURL key={key} action={value} data={data}/>
      }[value.type]
      if (component) {
        actionComponents.push(component)
      }
    }
    if (resource === "entries") {
      actionComponents.push(<EntryRowActions key="entries-default-action" data={data}/>)
    } else if (resource === "materials") {
      actionComponents.push(<MaterialRowActions key="materials-default-action" data={data}/>)
    }
    return actionComponents
  }, [rows, resource])

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
  let buttons
  if (resource === "entries") {
    details = rows?.details?.render || EntryDetails
    if (!noAction) buttons = <EntryDownloadButton tooltip="Download files" query={query} />
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
  columns: PropTypes.object,
  resource: PropTypes.string,
  rows: PropTypes.object,
  useResults: PropTypes.func,
  useApiQuery: PropTypes.func,
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

/**
 * Displays search results from the current search context.
 */
export const SearchResultsWithContext = React.memo((props) => {
  const {columns, resource, rows, useResults, useApiQuery} = useSearchContext()

  return <SearchResults
    columns={columns}
    resource={resource}
    rows={rows}
    useResults={useResults}
    useApiQuery={useApiQuery}
    {...props}
  />
})
