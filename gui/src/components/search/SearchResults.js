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
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper,
  Tooltip,
  IconButton
} from '@material-ui/core'
import ReloadIcon from '@material-ui/icons/Cached'
import { useQueryParam, useQueryParams, StringParam, NumberParam } from 'use-query-params'
import { searchContext } from './SearchContext'
import EntryList from './EntryList'
import ApiDialogButton from '../ApiDialogButton'
import {objectFilter} from '../../utils'

/**
 * Displays the list of search results
 */

const useStyles = makeStyles(theme => ({
  root: {}
}))

const SearchResults = React.memo(({
  availableTabs = ['entries', 'materials', 'datasets'],
  initialTab = 'entries',
  resultListProps = {},
  className
}) => {
  const styles = useStyles()
  const {domain} = useContext(searchContext)
  const ResultList = SearchEntryList

  return <div className={clsx(className, styles.root)}>
    <Paper>
      <ResultList domain={domain} {...resultListProps} />
    </Paper>
  </div>
})
SearchResults.propTypes = {
  availableTabs: PropTypes.arrayOf(PropTypes.string),
  initialTab: PropTypes.string,
  resultListProps: PropTypes.object,
  className: PropTypes.string
}

function ReRunSearchButton() {
  const {update} = useContext(searchContext)
  return <Tooltip title="Re-execute the search">
    <IconButton onClick={update}>
      <ReloadIcon />
    </IconButton>
  </Tooltip>
}

const usePagination = () => {
  const {setRequestParameters} = useContext(searchContext)
  let [requestQueryParameters, setRequestQueryParameters] = useQueryParams({
    order: NumberParam, order_by: StringParam, per_page: NumberParam, page: NumberParam
  })
  requestQueryParameters = objectFilter(requestQueryParameters, key => requestQueryParameters[key])
  requestQueryParameters.page = requestQueryParameters.page || 1
  useEffect(
    () => setRequestParameters(requestQueryParameters),
    [requestQueryParameters, setRequestParameters]
  )
  return setRequestQueryParameters
}

function SearchEntryList(props) {
  const {response, requestParameters, apiQuery, update} = useContext(searchContext)
  const setRequestParameters = usePagination()
  return <EntryList
    query={apiQuery}
    editable={apiQuery.owner === 'staging' || apiQuery.owner === 'user'}
    data={response}
    onChange={setRequestParameters}
    onEdit={update}
    actions={
      <>
        <ReRunSearchButton/>
        <ApiDialogButton data={response} />
      </>
    }
    {...requestParameters}
    {...props}
  />
}

export default SearchResults
