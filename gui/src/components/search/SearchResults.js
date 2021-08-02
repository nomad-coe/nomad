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
import React, { useCallback, useEffect, useState, useMemo, useRef } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper
} from '@material-ui/core'
import NewEntryList from './NewEntryList'
import { debounce } from 'lodash'
import { useApi } from '../apiV1'
import { useQuery, useExclusive, cleanQuery } from './FilterContext'

/**
 * Displays the list of search results
 */
const INIT_PAGE_SIZE = 30
const PAGE_SIZE_INCREMENT = 30

const useStyles = makeStyles(theme => ({
  root: {
    height: '100%'
  }
}))
const SearchResults = React.memo(({
  className
}) => {
  const styles = useStyles()
  const api = useApi()
  const query = useQuery()
  const exclusive = useExclusive()
  const [pagination, setPagination] = useState({
    page: 1,
    order: 'desc',
    order_by: 'upload_time'
  })
  const pageSize = useRef(INIT_PAGE_SIZE)
  const immediate = useRef(false)
  const [results, setResults] = useState()
  const total = results?.pagination.total || 0

  // API call for fetching results
  const loadResults = useCallback((query, pagination, exclusive) => {
    const search = {
      owner: 'all',
      query: cleanQuery(query, exclusive),
      pagination: pagination
    }
    api.queryEntry(search)
      .then(data => {
        setResults(data)
      })
  }, [api])

  // Debounced API call version
  const loadResultsDebounced = useCallback(debounce(loadResults, 400), [])

  // When bottom of results reached, increase the page size and ask request
  // immediate refresh.
  const handleBottom = useCallback(() => {
    immediate.current = true
    if (pageSize.current < total) {
      pageSize.current = Math.min(pageSize.current + PAGE_SIZE_INCREMENT, total)
      setPagination(old => {
        const newPagination = {...old, page_size: pageSize.current}
        return newPagination
      })
    }
  }, [total])

  // Handle sorting changes in the pagination: TODO: currently broken
  const handleChange = useCallback(({order, order_by}) => {
    immediate.current = true
    // setPagination(old => {
    //   const newPagination = {...old, order: order, order_by: order_by}
    //   return newPagination
    // })
  }, [])

  useEffect(() => {
    pageSize.current = INIT_PAGE_SIZE
    immediate.current = false
  }, [query, exclusive])

  // If the pagination changes, we load results immediately. Otherwise a
  // debounce is used.
  useEffect(() => {
    pagination.page_size = pageSize.current
    if (immediate.current) {
      loadResults(query, pagination, exclusive)
    } else {
      loadResultsDebounced(query, pagination, exclusive)
    }
  }, [query, pagination, exclusive, loadResults, loadResultsDebounced])

  // The rendered component is memoized in its entirety. This makes sure that we
  // re-render the results list only when the actual results have changed, and
  // not just when the search query changes. Has a significant effect on
  // performance.
  const comp = useMemo(() => {
    return <Paper className={clsx(className, styles.root)}>
      {results && <NewEntryList
        query={results.query}
        editable={results.query.owner === 'staging' || results.query.owner === 'user'}
        data={results}
        page={results?.pagination?.page}
        per_page={results?.pagination?.page_size}
        order={results?.pagination?.order}
        order_by={results?.pagination?.order_by}
        onChange={handleChange}
        onBottom={handleBottom}
      />}
    </Paper>
  }, [className, styles.root, results, handleChange, handleBottom])

  return comp
})
SearchResults.propTypes = {
  initialTab: PropTypes.string,
  resultListProps: PropTypes.object,
  className: PropTypes.string
}

export default SearchResults
