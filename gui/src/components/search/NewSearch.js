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
import React, {useState} from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import FiltersPanel from './FiltersPanel'
import FiltersActive from './FiltersActive'
import NewSearchBar from './NewSearchBar'
import SearchResults from './SearchResults'
import SearchContext from './SearchContext'

const useNewSearchStyles = makeStyles(theme => {
  const filterWidth = 26
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
      display: 'flex',
      flexDirection: 'column',
      paddingTop: theme.spacing(3),
      paddingBottom: theme.spacing(3),
      paddingLeft: theme.spacing(3),
      paddingRight: theme.spacing(4)
    },
    resultList: {
      overflowY: 'auto',
      flexGrow: 1
    },
    rightColumn: {
      //flex: `0 0 ${0.5 * filterWidth}rem`
    },
    spacer: {
      flexGrow: 1
    },
    bar: {
      display: 'flex',
      flexGrow: 0,
      marginBottom: theme.spacing(3)
    },
    searchBar: {
      flexGrow: 1,
      // flexBasis: '50%',
      //maxWidth: '35rem'
    },
    spacerBar: {
      flex: `0 0 ${theme.spacing(3)}px`
    },
    // filtersBar: {
    //   flexGrow: 1,
    //   flexBasis: '50%'
    // }
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
  const [resultType, setResultType] = useState('entries')
  const [searchType, setSearchType] = useState('nomad')
  return <SearchContext query={query} initialQuery={initialQuery}>
    <div className={styles.root} {...rest}>
      <div className={styles.leftColumn}>
        <FiltersPanel
          resultType={resultType}
          onResultTypeChange={value => setResultType(value)}
        />
      </div>
      <div className={styles.center}>
        <div className={styles.bar}>
          <NewSearchBar
            searchType={searchType}
            className={styles.searchBar}
            onSearchTypeChanged={(event) => setSearchType(event.target.value)}
          />
          {/* <div className={styles.spacerBar}></div>
          <FiltersActive className={styles.filtersBar}/> */}
        </div>
        <div className={styles.resultList}>
          <SearchResults/>
        </div>
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

export default NewSearch
