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
import React, { useState, useCallback, useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Divider,
  Paper,
  Menu,
  Typography
} from '@material-ui/core'
import ArrowForwardIcon from '@material-ui/icons/ArrowForward'
import ArrowBackIcon from '@material-ui/icons/ArrowBack'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import ClearIcon from '@material-ui/icons/Clear'
import RefreshIcon from '@material-ui/icons/Refresh'
import Scrollable from '../../visualization/Scrollable'
import FilterSummary from '../FilterSummary'
import FilterSettings from './FilterSettings'
import { Actions, ActionHeader, Action } from '../../Actions'
import { useSearchContext } from '../SearchContext'
import { filterGroups } from '../FilterRegistry'
import { MoreVert } from '@material-ui/icons'

// The menu animations use a transition on the 'transform' property. Notice that
// animating 'transform' instead of e.g. the 'left' property is much more
// performant. We also hint the browser that the transform property will be
// animated using the 'will-change' property: this will pre-optimize the element
// for animation when possible (the recommendation is to remove/add it when
// needed, but in this case we keep it on constantly).
//
// The menu widths are hardcoded. We need the widths to perform the close
// animation. Another option would be to use useLayoutEffect to determine the
// sizes of the components dynamically, but this seems to be quite a bit less
// responsive compared to hardcoding the values.

export const filterMenuContext = React.createContext()

const useFilterMenuStyles = makeStyles(theme => {
  const width = 22
  return {
    root: {
      boxSizing: 'border-box',
      display: 'flex',
      position: 'relative',
      flexDirection: 'column',
      width: `${width}rem`,
      height: '100%',
      '-webkit-transform': 'none',
      transform: 'none',
      transition: 'transform 250ms',
      willChange: 'transform'
    },
    collapsed: {
      '-webkit-transform': `translateX(-${width - 4}rem)`,
      transform: `translateX(-${width - 4}rem)`
    }
  }
})

export const FilterMenu = React.memo(({
  selected,
  onSelectedChange,
  open,
  onOpenChange,
  collapsed,
  onCollapsedChange,
  className,
  children
}) => {
  const styles = useFilterMenuStyles()
  const [size, setSize] = useState('medium')

  const handleChange = useCallback((newValue) => {
    if (newValue !== selected) {
      onOpenChange(true)
    } else {
      onOpenChange(old => !old)
    }
    onSelectedChange && onSelectedChange(newValue)
  }, [selected, onSelectedChange, onOpenChange])

  return <div className={clsx(className, styles.root, collapsed && styles.collapsed)}>
    <filterMenuContext.Provider value={{
      selected: selected,
      onChange: handleChange,
      open: open,
      onOpenChange: onOpenChange,
      size: size,
      onSizeChange: setSize,
      collapsed: collapsed,
      onCollapsedChange: onCollapsedChange
    }}>
      {children}
    </filterMenuContext.Provider>
  </div>
})
FilterMenu.propTypes = {
  selected: PropTypes.string,
  onSelectedChange: PropTypes.func,
  open: PropTypes.bool,
  onOpenChange: PropTypes.func,
  collapsed: PropTypes.bool,
  onCollapsedChange: PropTypes.func,
  className: PropTypes.string,
  children: PropTypes.node
}

const useFilterMenuItemsStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
  return {
    root: {
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'column',
      width: '100%',
      height: '100%',
      zIndex: 3
    },
    header: {
      paddingTop: theme.spacing(0.5),
      paddingBottom: theme.spacing(0.5),
      paddingLeft: padding,
      paddingRight: padding,
      overflow: 'visible'
    },
    headerText: {
      display: 'flex',
      alignItems: 'center'
    },
    headerTextVertical: {
      display: 'flex',
      alignItems: 'center',
      position: 'absolute',
      right: '0.5rem',
      top: '3.5rem',
      height: '4rem',
      transform: 'rotate(90deg)'
    },
    menu: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      left: 0,
      backgroundColor: theme.palette.background.paper,
      boxSizing: 'border-box'
    },
    menuBorder: {
      boxShadow: `1px 0px 0px 0px ${theme.palette.action.selected}`
    },
    padding: {
      paddingTop: `${theme.spacing(1.5)}px`,
      paddingBottom: `${theme.spacing(1.5)}px`
    },
    button: {
      marginRight: 0,
      '-webkit-transform': 'none',
      transform: 'none',
      transition: 'transform 250ms',
      willChange: 'transform'
    },
    hidden: {
      display: 'none'
    },
    overflow: {
      overflow: 'visible'
    }
  }
})

export const FilterMenuItems = React.memo(({
  className,
  children
}) => {
  const { useResetFilters, useRefresh } = useSearchContext()
  const styles = useFilterMenuItemsStyles()
  const { open, onOpenChange, collapsed, onCollapsedChange } = useContext(filterMenuContext)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const isSettingsOpen = Boolean(anchorEl)
  const resetFilters = useResetFilters()
  const refresh = useRefresh()

  // Callbacks
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

  // Unfortunately the ClickAwayListener does not play nicely together with
  // Menus/Select/Popper. When using Portals, the clicks are registered wrong.
  // When Portals are disabled (disablePortal), their positioning goes haywire.
  // The clicks outside are thus detected by individual event listeners that
  // toggle the menu state.
  return <div className={clsx(className, styles.root)}>
    <div
      className={clsx(styles.menu, open && styles.menuBorder, collapsed && styles.hidden)}
      classes={{containerInner: styles.overflow}}>
      <Scrollable>
        <div className={styles.padding}>
          <Actions className={styles.header}>
            <ActionHeader>
              <Typography className={styles.headerText} variant="button">Filters</Typography>
            </ActionHeader>
            <Action
              tooltip="Refresh results"
              onClick={() => refresh()}
            >
              <RefreshIcon/>
            </Action>
            <Action
              tooltip="Clear filters"
              onClick={() => resetFilters()}
            >
              <ClearIcon/>
            </Action>
            <Action
              tooltip="Options"
              onClick={openMenu}
            >
              <MoreVert/>
            </Action>
            <Menu
              anchorEl={anchorEl}
              open={isSettingsOpen}
              onClose={closeMenu}
              getContentAnchorEl={null}
              anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
              transformOrigin={{ vertical: 'top', horizontal: 'right' }}
              keepMounted
            >
              <div>
                <FilterSettings/>
              </div>
            </Menu>
            {!collapsed && <Action
              tooltip={'Hide filter menu'}
              onClick={() => {
                onCollapsedChange(old => !old)
                onOpenChange(false)
              }}
              className={styles.button}
            >
              <ArrowBackIcon/>
            </Action>}
          </Actions>
          <List dense>
            <Divider/>
            {children}
          </List>
        </div>
      </Scrollable>
    </div>
    <div className={clsx(styles.padding, !collapsed && styles.hidden)}>
      <Actions className={styles.header}>
        {collapsed && <Action
          tooltip={'Show filter menu'}
          onClick={() => {
            onCollapsedChange(false)
          }}
          className={styles.button}
        >
          <ArrowForwardIcon/>
        </Action>}
      </Actions>
      <Typography
        className={clsx(styles.headerText, styles.headerTextVertical)}
        variant="button"
      >Filters
      </Typography>
    </div>
  </div>
})
FilterMenuItems.propTypes = {
  className: PropTypes.string,
  children: PropTypes.node
}

const useFilterMenuItemStyles = makeStyles(theme => {
  return {
    root: {},
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
    listItem: {
      position: 'relative',
      height: '2.6rem'
    },
    divider: {
      width: '100%',
      backgroundColor: theme.palette.grey[300]
    }
  }
})
export const FilterMenuItem = React.memo(({
  value,
  group,
  onClick,
  disableButton,
  depth
}) => {
  const styles = useFilterMenuItemStyles()
  const theme = useTheme()
  const groupFinal = group || filterGroups[value]
  const { selected, open, onChange } = useContext(filterMenuContext)
  const handleClick = disableButton ? undefined : (onClick || onChange)
  const opened = open && value === selected

  return <>
    <ListItem
      button={!!handleClick}
      className={styles.listItem}
      classes={{gutters: styles.gutters}}
      onClick={handleClick && (() => handleClick(value))}
    >
      <ListItemText
        style={{marginLeft: theme.spacing(depth * 1.8)}}
        primaryTypographyProps={{
          color: opened ? 'primary' : 'initial',
          className: styles.label
        }}
        primary={value}
      />
      {handleClick && <ListItemIcon className={styles.listIcon}>
        <NavigateNextIcon color={opened ? 'primary' : 'action'} className={styles.arrow}/>
      </ListItemIcon>}
    </ListItem>
    {groupFinal && <FilterSummary quantities={groupFinal}/>}
    <Divider className={styles.divider}/>
  </>
})

FilterMenuItem.propTypes = {
  value: PropTypes.string,
  group: PropTypes.string,
  onClick: PropTypes.func,
  disableButton: PropTypes.bool,
  depth: PropTypes.number
}
FilterMenuItem.defaultProps = {
  depth: 0
}

const useFilterSubMenuStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  },
  hidden: {
    display: 'none'
  }
}))

const useFilterSubMenusStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
  const widthMedium = 25
  const widthLarge = 48
  return {
    root: {
      zIndex: 2,
      display: 'flex',
      flexDirection: 'column',
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      width: `${widthLarge}rem`,
      backgroundColor: theme.palette.background.paper,
      '-webkit-transform': 'none',
      transform: 'none',
      transition: 'transform 250ms',
      flexGrow: 1,
      boxSizing: 'border-box',
      willChange: 'transform'
    },
    collapsed: {
      display: 'none'
    },
    header: {
      paddingTop: theme.spacing(0.5),
      paddingBottom: theme.spacing(1.5),
      paddingRight: theme.spacing(0),
      paddingLeft: theme.spacing(0)
    },
    headerText: {
      display: 'flex',
      alignItems: 'center'
    },
    containerMedium: {
      '-webkit-transform': `translateX(${widthMedium}rem)`,
      transform: `translateX(${widthMedium}rem)`
    },
    containerLarge: {
      '-webkit-transform': `translateX(${widthLarge}rem)`,
      transform: `translateX(${widthLarge}rem)`
    },
    menu: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      boxSizing: 'border-box'
    },
    padding: {
      padding: `${theme.spacing(1.5)}px ${padding}px`
    },
    menuMedium: {
      width: `${widthMedium}rem`
    },
    menuLarge: {
      width: `${widthLarge}rem`
    },
    button: {
      marginRight: 0
    }
  }
})

export const FilterSubMenus = React.memo(({
  children
}) => {
  const styles = useFilterSubMenusStyles()
  const { selected, open, onOpenChange, size, collapsed } = useContext(filterMenuContext)
  const [menuStyle, containerStyle] = {
    medium: [styles.menuMedium, styles.containerMedium],
    large: [styles.menuLarge, styles.containerLarge]
  }[size]

  return <Paper
    elevation={4}
    className={clsx(styles.root, open && containerStyle)}
  >
    <div className={clsx(styles.menu, menuStyle, collapsed && styles.collapsed)}>
      <Scrollable>
        <div className={styles.padding}>
          <Actions className={styles.header}>
            <ActionHeader>
              <Typography className={styles.headerText} variant="button">{selected}</Typography>
            </ActionHeader>
            <Action
              tooltip="Close submenu"
              onClick={() => { onOpenChange(false) }}
              className={styles.button}
            >
              <ArrowBackIcon/>
            </Action>
          </Actions>
          {children}
        </div>
      </Scrollable>
    </div>
  </Paper>
})
FilterSubMenus.propTypes = {
  sizes: PropTypes.object,
  children: PropTypes.node
}

export const FilterSubMenu = React.memo(({
  value,
  size,
  children
}) => {
  const styles = useFilterSubMenuStyles()
  const { selected, onSizeChange } = useContext(filterMenuContext)
  const visible = value === selected
  useEffect(() => {
    if (visible) {
      onSizeChange(size)
    }
  }, [size, visible, onSizeChange])

  return <div className={clsx(styles.root, !visible && styles.hidden)}>
    {children}
  </div>
})
FilterSubMenu.propTypes = {
  value: PropTypes.string,
  size: PropTypes.oneOf(['medium', 'large']),
  children: PropTypes.node
}
FilterSubMenu.defaultProps = {
  size: 'medium'
}

export default FilterSubMenu
