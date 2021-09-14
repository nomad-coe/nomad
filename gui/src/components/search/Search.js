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
import React, { useState } from 'react'
import clsx from 'clsx'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import FilterMainMenu from './menus/FilterMainMenu'
import SearchBar from './SearchBar'
import SearchResults from './results/SearchResults'
import {
  useMenuOpenState
} from './SearchContext'

/**
 * The primary search interface that is reused throughout the application in
 * different contexts. Displays a menu of filters, a search bar, a list of
 * results and optionally a customizable header above the search bar.
 */
const useStyles = makeStyles(theme => {
  return {
    root: {
      display: 'flex',
      height: '100%',
      width: '100%',
      overflow: 'hidden'
    },
    leftColumn: {
      flexShrink: 0,
      flexGrow: 0,
      height: '100%',
      zIndex: 2
    },
    leftColumnCollapsed: {
      maxWidth: '4rem'
    },
    center: {
      flex: `1 1 100%`,
      display: 'flex',
      flexDirection: 'column',
      zIndex: 1,
      paddingBottom: theme.spacing(2.5),
      paddingLeft: theme.spacing(3),
      paddingRight: theme.spacing(3),
      paddingTop: theme.spacing(0.25)
    },
    container: {
    },
    resultList: {
      flexGrow: 1,
      minHeight: 0 // This makes sure that the flex item is not bigger than the parent.
    },
    spacer: {
      flexGrow: 1
    },
    header: {
      marginTop: theme.spacing(2)
    },
    searchBar: {
      marginTop: theme.spacing(2),
      display: 'flex',
      flexGrow: 0,
      zIndex: 1,
      marginBottom: theme.spacing(2.0)
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
      willChange: 'opacity',
      zIndex: 1
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

const Search = React.memo(({
  collapsed,
  header
}) => {
  const styles = useStyles()
  const [isMenuOpen, setIsMenuOpen] = useMenuOpenState(false)
  const [isCollapsed, setIsCollapsed] = useState(collapsed)

  return <div className={styles.root}>
    <div className={clsx(styles.leftColumn, isCollapsed && styles.leftColumnCollapsed)}>
      <FilterMainMenu
        open={isMenuOpen}
        onOpenChange={setIsMenuOpen}
        collapsed={isCollapsed}
        onCollapsedChange={setIsCollapsed}
      />
    </div>
    <div className={styles.center} onClick={() => setIsMenuOpen(false)}>
      <div className={styles.header}>
        {header}
      </div>
      <SearchBar
        className={styles.searchBar}
      />
      <SearchResults
        className={styles.resultList}
      />
      <div className={clsx(styles.nonInteractive, styles.shadow, styles.shadowHidden, isMenuOpen && styles.shadowVisible)}></div>
    </div>
  </div>
})
Search.propTypes = {
  collapsed: PropTypes.bool,
  header: PropTypes.node
}

export default Search
