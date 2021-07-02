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
  return {
    root: {},
    list: {
      paddingTop: 0,
      paddingBottom: 0
    },
    listChild: {
      marginLeft: theme.spacing(1.8)
    },
    label: {
      textTransform: 'capitalize'
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
      paddingLeft: theme.spacing(3.0),
      paddingRight: theme.spacing(2.35)
    },
    li: {
      display: 'flex',
      flexDirection: 'column',
      width: '100%'
    },
    listItem: {
      position: 'relative',
      height: '2.6rem'
    },
    divider: {
      width: '100%'
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
      const component = item.component
      const isOpen = isMenuOpen && view === childName

      // Add this item
      list.push(<div
        key={list.length}
        className={styles.li}
      >
        <ListItem
          button={!!component}
          className={styles.listItem}
          classes={{gutters: styles.gutters}}
          onClick={component ? () => handleClick(childName) : undefined}
        >
          <ListItemText
            className={level === 0 ? undefined : styles.listChild}
            primaryTypographyProps={{
              color: isOpen ? 'primary' : 'initial',
              className: styles.label
            }}
            primary={childName}
          />
          {component && <ListItemIcon className={styles.listIcon}>
            <NavigateNextIcon color={isOpen ? 'primary' : 'action'} className={styles.arrow}/>
          </ListItemIcon>}
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
