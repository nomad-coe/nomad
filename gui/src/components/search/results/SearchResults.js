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
import React, { useCallback, useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper
} from '@material-ui/core'
import SearchResultsMaterials from './SearchResultsMaterials'
import SearchResultsEntries from './SearchResultsEntries'
import { useExclusive, useScrollResults, useSearchContext } from '../FilterContext'

/**
 * Displays the list of search results
 */
const useStyles = makeStyles(theme => ({
  root: {
    height: '100%'
  }
}))
const orderByMap = {
  'entries': 'upload_time',
  'materials': 'chemical_formula_hill'
}
const SearchResults = React.memo(({
  className
}) => {
  const styles = useStyles()
  const exclusive = useExclusive()
  const { resource } = useSearchContext()
  const page_size = 30
  // eslint-disable-next-line no-unused-vars
  const [orderBy, setOrderBy] = useState(orderByMap[resource])
  // eslint-disable-next-line no-unused-vars
  const [order, setOrder] = useState('desc')
  const {results, next, page, total} = useScrollResults(30, orderBy, order, exclusive)

  // When bottom of results reached, fetch the next set of results.
  const handleBottom = useCallback(() => {
    if (results.data.length < total) {
      next()
    }
  }, [results, total, next])

  // Handle sorting changes in the pagination: TODO: currently broken.
  const handleChange = useCallback(({order, order_by}) => {
  }, [])

  // The rendered component is memoized in its entirety. This makes sure that we
  // re-render the results list only when the actual results have changed, and
  // not just when the search query changes. Has a significant effect on
  // performance.
  // const component = resource === 'materials' ? MaterialResults : NewEntryList
  const result = useMemo(() => {
    const Component = resource === 'materials' ? SearchResultsMaterials : SearchResultsEntries
    return <Paper className={clsx(className, styles.root)}>
      {results && <Component
        query={results.query}
        editable={false}
        data={results}
        page={page}
        per_page={page_size}
        order={order}
        order_by={orderBy}
        onChange={handleChange}
        onBottom={handleBottom}
      />}
    </Paper>
  }, [resource, className, styles.root, results, page, order, orderBy, handleChange, handleBottom])

  return result
})
SearchResults.propTypes = {
  className: PropTypes.string
}

export default SearchResults
