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
import React, { useEffect, useState } from 'react'
import clsx from 'clsx'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import FilterMainMenu from './menus/FilterMainMenu'
import NewSearchBar from './NewSearchBar'
import SearchResults from './results/SearchResults'
import {
  quantities,
  useMenuOpenState,
  useInitQuery,
  useInitialAggs,
  useSetOwner
} from './FilterContext'

const useStyles = makeStyles(theme => {
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
      flexGrow: 1,
      minHeight: 0 // This makes sure that the flex item is not bigger than the parent.
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
      zIndex: 1
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

const NewSearch = React.memo(({
  owner,
  collapsed
}) => {
  const styles = useStyles()
  const [resultType, setResultType] = useState('entries')
  const [isMenuOpen, setIsMenuOpen] = useMenuOpenState(false)
  const [isCollapsed, setIsCollapsed] = useState(collapsed)
  const setOwner = useSetOwner()
  useInitQuery()
  useInitialAggs()

  useEffect(() => {
    setOwner(owner)
  }, [setOwner, owner])

  return <div className={styles.root}>
    <div className={styles.leftColumn}>
      <FilterMainMenu
        open={isMenuOpen}
        resultType={resultType}
        onResultTypeChange={value => setResultType(value)}
        onOpenChange={setIsMenuOpen}
        collapsed={isCollapsed}
        onCollapseChanged={setIsCollapsed}
      />
    </div>
    <div className={styles.center} onClick={ () => { setIsMenuOpen(false) } }>
      <div className={styles.bar}>
        <NewSearchBar
          quantities={quantities}
          className={styles.searchBar}
        />
      </div>
      <div className={styles.resultList}>
        <SearchResults/>
      </div>
      <div className={clsx(styles.nonInteractive, styles.shadow, styles.shadowHidden, isMenuOpen && styles.shadowVisible)}></div>
    </div>
  </div>
})
NewSearch.propTypes = {
  owner: PropTypes.string,
  collapsed: PropTypes.bool
}

export default NewSearch
