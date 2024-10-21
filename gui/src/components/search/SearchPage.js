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
import React from 'react'
import clsx from 'clsx'
import PropTypes from 'prop-types'
import { Box } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import SearchMenu from './SearchMenu'
import SearchBar from './SearchBar'
import Query from './Query.js'
import { SearchResults } from './SearchResults'
import Dashboard from './widgets/Dashboard'
import { useSearchContext } from './SearchContext'

/**
 * The primary search interface that is reused throughout the application in
 * different contexts. Displays a menu of filters, a search bar, a list of
 * results and optionally a customizable header above the search bar.
 */
const useStyles = makeStyles(theme => {
  return {
    root: {
      display: 'flex',
      flexDirection: 'row',
      height: '100%',
      width: '100%'
    },
    leftColumn: {
      flexShrink: 0,
      flexGrow: 0,
      height: '100%',
      zIndex: 2
    },
    center: {
      flexGrow: 1,
      height: '100%',
      overflowY: 'scroll'
    },
    searchBar: {
      display: 'flex',
      flexGrow: 0,
      zIndex: 1,
      marginBottom: theme.spacing(0.3)
    },
    shadow: {
      pointerEvents: 'none',
      position: 'absolute',
      left: 0,
      right: 0,
      bottom: 0,
      top: 0,
      height: '100%',
      width: '100%',
      backgroundColor: 'black',
      transition: 'opacity 200ms',
      willChange: 'opacity',
      zIndex: 1,
      opacity: 0
    },
    shadowVisible: {
      opacity: 0.1
    }
  }
})

const SearchPage = React.memo(({
  header
}) => {
  const styles = useStyles()
  const {useIsMenuOpen, useSetIsMenuOpen} = useSearchContext()
  const [isMenuOpen, setIsMenuOpen] = [useIsMenuOpen(), useSetIsMenuOpen()]

  return <div className={styles.root}>
    <div className={clsx(styles.leftColumn)}>
      <SearchMenu/>
    </div>
    <div className={styles.center} onClick={() => setIsMenuOpen(false)}>
      <Box margin={2.5} paddingBottom={3}>
        <Box marginBottom={2}>
          {header}
        </Box>
        <Box marginBottom={0}>
          <SearchBar className={styles.searchBar} />
          <Query/>
        </Box>
        <Box marginBottom={1} zIndex={0}>
          <Dashboard/>
        </Box>
        <Box position="relative" zIndex={1}>
          <SearchResults/>
        </Box>
        <div className={clsx(styles.shadow, isMenuOpen && styles.shadowVisible)}></div>
      </Box>
    </div>
  </div>
})
SearchPage.propTypes = {
  header: PropTypes.node
}

export default SearchPage
