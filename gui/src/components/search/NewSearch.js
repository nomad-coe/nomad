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
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper,
  Tabs,
  Tab,
  Tooltip,
  IconButton
} from '@material-ui/core'
import ReloadIcon from '@material-ui/icons/Cached'
import { useQueryParam, useQueryParams, StringParam, NumberParam } from 'use-query-params'
import SearchBar from './SearchBar'
import FiltersPanel from './FiltersPanel'
import EntryList from './EntryList'
import DatasetList from './DatasetList'
import MaterialsList from './MaterialsList'
import ApiDialogButton from '../ApiDialogButton'
import SearchContext, { searchContext } from './SearchContext'
import {objectFilter} from '../../utils'

const resultTabs = {
  'entries': {
    label: 'Entries',
    groups: {},
    component: SearchEntryList
  },
  'materials': {
    label: 'Materials',
    groups: {'encyclopedia.material.materials_grouped': true},
    component: SearchMaterialsList
  },
  'datasets': {
    label: 'Datasets',
    groups: {'datasets_grouped': true},
    component: SearchDatasetList
  }
}

const useNewSearchStyles = makeStyles(theme => {
  const filterWidth = 25
  return {
    root: {
      display: 'flex',
      height: '100%',
      width: '100%',
      overflow: 'hidden'
    },
    leftColumn: {
      flex: `0 0 ${filterWidth}rem`,
      maxWidth: `${filterWidth}rem`,
      height: '100%',
      position: 'relative'
    },
    center: {
      flex: `1 1 100%`,
      padding: theme.spacing(3)
    },
    rightColumn: {
      flex: `0 0 ${0.5 * filterWidth}rem`
    },
    spacer: {
      flexGrow: 1
    },
    searchBar: {
      marginTop: theme.spacing(1),
      marginBottom: theme.spacing(1)
    }
  }
})

const NewSearch = React.memo(({
  initialOwner,
  ownerTypes,
  initialMetric,
  query,
  initialQuery,
  resultListProps,
  initialRequest,
  showDisclaimer,
  ...rest
}) => {
  const styles = useNewSearchStyles()
  return <SearchContext query={query} initialQuery={initialQuery}>
    <div className={styles.root} {...rest}>
      <div className={styles.leftColumn}>
        <FiltersPanel/>
      </div>
      <div className={styles.center}>
        <SearchBar classes={{autosuggestRoot: styles.searchBar}}/>
        <SearchResults/>
      </div>
      <div className={styles.rightColumn}>
      </div>
    </div>
  </SearchContext>
})
NewSearch.propTypes = {
  initialOwner: PropTypes.string,
  ownerTypes: PropTypes.arrayOf(PropTypes.string),
  initialMetric: PropTypes.string,
  initialRequest: PropTypes.object,
  resultListProps: PropTypes.object,
  /**
   * Additional search parameters that will be added to all searches that are send to
   * the API. The idea is that this can be used to lock some aspects of the search for
   * special contexts, like the dataset page for example.
   */
  query: PropTypes.object,
  /**
   * Similar to query, but these parameters can be changes by the user interacting with
   * the component.
   */
  initialQuery: PropTypes.object,
  showDisclaimer: PropTypes.bool
}

const useSearchResultStyles = makeStyles(theme => ({
  root: {}
}))
const SearchResults = React.memo(({
  availableTabs = ['entries'],
  initialTab = 'entries',
  resultListProps = {}
}) => {
  const classes = useSearchResultStyles()
  const {domain, setGroups} = useContext(searchContext)
  let [openTab, setOpenTab] = useQueryParam('results', StringParam)
  openTab = openTab || initialTab
  const ResultList = resultTabs[openTab].component
  const handleTabChange = tab => {
    setOpenTab(tab)
    setGroups(resultTabs[tab].groups)
  }

  useEffect(() => {
    if (openTab !== 'entries') {
      handleTabChange(openTab)
    }
    // eslint-disable-next-line
  }, [])

  return <div className={classes.root}>
    <Paper>
      <Tabs
        value={openTab}
        indicatorColor="primary"
        textColor="primary"
        onChange={(event, value) => handleTabChange(value)}
      >
        {availableTabs.filter(tab => domain.searchTabs.includes(tab)).map(key => {
          const tab = resultTabs[key]
          return <Tab key={key} label={tab.label} value={key} />
        })}
      </Tabs>

      <ResultList domain={domain} {...resultListProps} />
    </Paper>
  </div>
})
SearchResults.propTypes = {
  'availableTabs': PropTypes.arrayOf(PropTypes.string),
  'initialTab': PropTypes.string,
  'resultListProps': PropTypes.object
}

export default NewSearch

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

const useScroll = (apiGroupName, afterParameterName) => {
  afterParameterName = afterParameterName || `${apiGroupName}_after`
  const apiAfterParameterName = `${apiGroupName}_grouped_after`

  const {response, setRequestParameters} = useContext(searchContext)
  const [queryAfterParameter, setQueryAfterParameter] = useQueryParam(afterParameterName, StringParam)
  useEffect(
    () => {
      const requestParameters = {}
      requestParameters[apiAfterParameterName] = queryAfterParameter || null
      setRequestParameters(requestParameters)
    }, [queryAfterParameter, setRequestParameters, apiAfterParameterName]
  )

  const responseGroup = response[`${apiGroupName}_grouped`]
  const after = responseGroup && responseGroup.after
  const result = {
    total: response.statistics.total.all[apiGroupName],
    onChange: requestParameters => setQueryAfterParameter(requestParameters[apiAfterParameterName])
  }
  result[afterParameterName] = after
  return result
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
      <React.Fragment>
        <ReRunSearchButton/>
        <ApiDialogButton data={response} />
      </React.Fragment>
    }
    {...requestParameters}
    {...props}
  />
}

function SearchDatasetList(props) {
  const {response, update} = useContext(searchContext)
  return <DatasetList
    data={response}
    onEdit={update}
    actions={<ReRunSearchButton/>}
    {...response} {...props} {...useScroll('datasets')}
  />
}

function SearchMaterialsList(props) {
  const {response} = useContext(searchContext)
  return <MaterialsList
    data={response}
    actions={<ReRunSearchButton/>}
    {...response} {...props} {...useScroll('encyclopedia.material.materials', 'materials_after')}
  />
}
