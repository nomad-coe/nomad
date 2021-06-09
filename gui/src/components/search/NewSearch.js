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
import clsx from 'clsx'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import FilterPanel from './FilterPanel'
import NewSearchBar from './NewSearchBar'
import SearchResults from './SearchResults'
import SearchContext from './SearchContext'
import { searchBarQuantities } from './FilterContext'

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
      position: 'relative',
      flex: `1 1 100%`,
      display: 'flex',
      flexDirection: 'column',
      paddingTop: theme.spacing(3),
      paddingBottom: theme.spacing(3),
      paddingLeft: theme.spacing(3),
      paddingRight: theme.spacing(4)
    },
    resultList: {
      // overflowY: 'auto',
      flexGrow: 1
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
      flexGrow: 1
    },
    spacerBar: {
      flex: `0 0 ${theme.spacing(3)}px`
    },
    nonInteractive: {
      pointerEvents: 'none',
      position: 'absolute',
      left: 0,
      right: 0,
      bottom: 0,
      top: 0,
      height: '100%',
      width: '100%'
    },
    shadow: {
      backgroundColor: 'black',
      transition: 'opacity 200ms',
      willChange: 'opacity'
    },
    hidden: {
      display: 'none'
    },
    shadowHidden: {
      opacity: 0
    },
    shadowVisible: {
      opacity: 0.1
    },
    placeholderVisible: {
      display: 'block'
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
  const [resultType, setResultType] = useState('entries')
  const [isMenuOpen, setIsMenuOpen] = useState(false)
  return <SearchContext query={query} initialQuery={initialQuery}>
    <div className={styles.root} {...rest}>
      <div className={styles.leftColumn}>
        <FilterPanel
          isMenuOpen={isMenuOpen}
          resultType={resultType}
          onResultTypeChange={value => setResultType(value)}
          onIsMenuOpenChange={setIsMenuOpen}
        />
      </div>
      <div className={styles.center}>
        <div className={styles.bar}>
          <NewSearchBar
            quantities={searchBarQuantities}
            className={styles.searchBar}
          />
        </div>
        <div className={styles.resultList}>
          <SearchResults/>
        </div>
        <div className={clsx(styles.nonInteractive, styles.shadow, styles.shadowHidden, isMenuOpen && styles.shadowVisible)}></div>
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
