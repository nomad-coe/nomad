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
import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper
} from '@material-ui/core'
import NewEntryList from './NewEntryList'
import { useResults } from './FilterContext'

/**
 * Displays the list of search results
 */
const useStyles = makeStyles(theme => ({
  root: {
    height: '100%'
  }
}))
const SearchResults = React.memo(({
  className
}) => {
  const styles = useStyles()

  // Find out how many results will fit on the current page.
  const [pagination, setPagination] = useState({
    page: 1,
    page_size: 20,
    order: 'desc',
    order_by: 'upload_time'
  })

  const results = useResults(pagination)
  const total = results?.pagination.total || 0

  // Request new results when scrolled to bottom. TODO: we should definitely use
  // scrolling here instead of increasing the page size...
  const handleBottom = useCallback(() => {
    setPagination(old => {
      if (old.page_size >= total) {
        return old
      }
      const newPagination = {...old, page_size: Math.min(old.page_size + 20, total)}
      return newPagination
    })
  }, [total])

  return <Paper className={clsx(className, styles.root)}>
    {results && <NewEntryList
      query={results.query}
      editable={results.query.owner === 'staging' || results.query.owner === 'user'}
      data={results}
      page={results?.pagination?.page}
      per_page={results?.pagination?.page_size}
      order={results?.pagination?.order}
      order_by={results?.pagination?.order_by}
      onBottom={handleBottom}
    />}
  </Paper>
})
SearchResults.propTypes = {
  initialTab: PropTypes.string,
  resultListProps: PropTypes.object,
  className: PropTypes.string
}

export default SearchResults
