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
import React, { useMemo, useCallback } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import {
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Divider
} from '@material-ui/core'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'

/**
 * Displays the tree-like structure for selecting filters.
 */
const useStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
  return {
    root: {},
    list: {
      paddingTop: 0,
      paddingBottom: 0
    },
    section: {
    },
    hidden: {
      display: 'none'
    },
    listChild: {
      marginLeft: theme.spacing(1)
    },
    listParent: {
    },
    listIcon: {
      fontsize: '1rem',
      minWidth: '1.5rem'
    },
    arrow: {
      marginLeft: theme.spacing(1),
      fontSize: '1.5rem'
    },
    gutters: {
      paddingLeft: theme.spacing(3),
      paddingRight: padding
    },
    li: {
      display: 'flex',
      flexDirection: 'column',
      width: '100%'
    },
    listItem: {
      height: '2.5rem'
    },
    divider: {
      width: '100%'
    },
    selected: {
      '&$selected': {
        backgroundColor: theme.palette.primary.main,
        color: 'black'
      },
      '&$selected:hover': {
        backgroundColor: theme.palette.primary.main,
        color: 'black'
      }
    }
  }
})

const FiltersTree = React.memo(({
  filterTree,
  view,
  isMenuOpen,
  onViewChange,
  onIsMenuOpenChange,
  className
}) => {
  const styles = useStyles()

  // Handle menu item click
  const handleClick = useCallback(view => {
    onViewChange(oldView => {
      if (view !== oldView) {
        onIsMenuOpenChange(true)
      } else {
        onIsMenuOpenChange(oldOpen => !oldOpen)
      }
      return view
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Determine the navigation tree layout
  const tree = useMemo(() => {
    function buildList(list, item, level) {
      const childName = item.name
      const filters = item.filters
      const isOpen = isMenuOpen && view === childName

      // Add this item
      list.push(<div
        key={list.length}
        className={styles.li}
      >
        <ListItem
          button
          className={styles.listItem}
          classes={{gutters: styles.gutters}}
          onClick={() => handleClick(childName)}
        >
          <ListItemText
            className={level === 0 ? styles.listParent : styles.listChild}
            primaryTypographyProps={{color: isOpen ? 'primary' : 'textPrimary'}}
            primary={childName}
          />
          <ListItemIcon className={styles.listIcon}>
            <NavigateNextIcon color={isOpen ? 'primary' : 'action'} className={styles.arrow}/>
          </ListItemIcon>
        </ListItem>
        {filters}
        <Divider className={styles.divider}/>
      </div>)

      // Add children recursively
      const children = item.children
      if (children) {
        for (let child of children) {
          buildList(list, child, level + 1)
        }
      }
    }
    const itemList = []
    for (let item of filterTree) {
      buildList(itemList, item, 0)
    }
    return <List dense className={styles.list}>
      <Divider/>
      {itemList}
    </List>
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [view, isMenuOpen])

  return <div className={clsx(className, styles.root)}>
    {tree}
  </div>
})
FiltersTree.propTypes = {
  filterTree: PropTypes.array,
  view: PropTypes.string,
  isMenuOpen: PropTypes.bool,
  className: PropTypes.string,
  onViewChange: PropTypes.func,
  onIsMenuOpenChange: PropTypes.func
}
FiltersTree.defaultProps = {
  level: 0
}

export default FiltersTree
