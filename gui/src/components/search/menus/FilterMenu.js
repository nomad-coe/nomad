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
import MoreVert from '@material-ui/icons/MoreVert'
import Scrollable from '../../visualization/Scrollable'
import FilterSummary from '../FilterSummary'
import FilterSettings from './FilterSettings'
import { Actions, ActionHeader, Action } from '../../Actions'
import { useSearchContext } from '../SearchContext'
import { filterGroups } from '../FilterRegistry'
import { pluralize } from '../../../utils'
import { isNil } from 'lodash'

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
const paddingHorizontal = 1.5
export const collapsedMenuWidth = 3.3

// Topmost header for filter menus. Contains menu actions and an overline text
const useFilterMenuTopHeaderStyles = makeStyles(theme => {
  return {
    root: {
      paddingTop: theme.spacing(0.7),
      paddingRight: theme.spacing(paddingHorizontal),
      paddingLeft: theme.spacing(paddingHorizontal)
    },
    title: {
      display: 'flex',
      alignItems: 'center',
      fontSize: '0.75rem',
      marginBottom: -2
    }
  }
})
export const FilterMenuTopHeader = React.memo(({
  title,
  actions,
  className
}) => {
  const styles = useFilterMenuTopHeaderStyles()
  return <div className={clsx(className, styles.root)}>
    <Actions>
      <ActionHeader>
        <Typography className={styles.title} variant="overline">{title}</Typography>
      </ActionHeader>
      {actions}
    </Actions>
  </div>
})
FilterMenuTopHeader.propTypes = {
  overlineTitle: PropTypes.string,
  title: PropTypes.string,
  topAction: PropTypes.node,
  actions: PropTypes.node,
  className: PropTypes.string
}

// Header for filter menus. Contains panel actions and an overline text.
const useFilterMenuHeaderStyles = makeStyles(theme => {
  return {
    root: {
      height: theme.spacing(3),
      paddingBottom: theme.spacing(0.8),
      paddingTop: theme.spacing(0.30),
      paddingLeft: theme.spacing(paddingHorizontal),
      paddingRight: theme.spacing(paddingHorizontal),
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'center'
    },
    title: {
      display: 'flex',
      alignItems: 'center',
      fontSize: '0.90rem'
    }
  }
})
export const FilterMenuHeader = React.memo(({
  title,
  actions,
  className
}) => {
  const styles = useFilterMenuHeaderStyles()
  return <div className={clsx(className, styles.root)}>
    <Actions>
      <ActionHeader>
        <Typography className={styles.title} variant="button">{title}</Typography>
      </ActionHeader>
      {actions}
    </Actions>
  </div>
})

FilterMenuHeader.propTypes = {
  overlineTitle: PropTypes.string,
  title: PropTypes.string,
  topAction: PropTypes.node,
  actions: PropTypes.node,
  className: PropTypes.string
}

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
      '-webkit-transform': `translateX(-${width - collapsedMenuWidth}rem)`,
      transform: `translateX(-${width - collapsedMenuWidth}rem)`
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
  return {
    root: {
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'column',
      width: '100%',
      height: '100%',
      backgroundColor: theme.palette.background.paper,
      zIndex: 3
    },
    headerTextVertical: {
      display: 'flex',
      alignItems: 'center',
      position: 'absolute',
      right: '0.07rem',
      top: '2.5rem',
      height: `${collapsedMenuWidth}rem`,
      transform: 'rotate(90deg)'
    },
    menu: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      left: 0,
      boxSizing: 'border-box'
    },
    menuBorder: {
      boxShadow: `1px 0px 0px 0px ${theme.palette.action.selected}`
    },
    padding: {
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
    },
    list: {
      paddingTop: 0
    },
    container: {
      display: 'flex',
      flexDirection: 'column',
      height: '100%'
    },
    content: {
      flex: 1,
      minHeight: 0
    },
    collapsedActions: {
      paddingTop: theme.spacing(0.7),
      paddingLeft: theme.spacing(paddingHorizontal),
      paddingRight: theme.spacing(paddingHorizontal)
    }
  }
})
export const FilterMenuItems = React.memo(({
  className,
  children
}) => {
  const { useResetFilters, useRefresh, resource } = useSearchContext()
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
    <div className={clsx(styles.menu, open && styles.menuBorder, collapsed && styles.hidden)}>
      <div className={styles.container}>
        <FilterMenuTopHeader
          title={`${resource} search`}
          actions={!collapsed && <Action
            tooltip={'Hide filter menu'}
            onClick={() => {
              onCollapsedChange(old => !old)
              onOpenChange(false)
            }}
            className={styles.button}
          >
            <ArrowBackIcon fontSize="small"/>
          </Action>}
        />
        <FilterMenuHeader
          title="Filters"
          actions={<>
            <Action
              tooltip="Refresh results"
              onClick={() => refresh()}
            >
              <RefreshIcon fontSize="small"/>
            </Action>
            <Action
              tooltip="Clear filters"
              onClick={() => resetFilters()}
            >
              <ClearIcon fontSize="small"/>
            </Action>
            <Action
              tooltip="Options"
              onClick={openMenu}
            >
              <MoreVert fontSize="small"/>
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
          </>}
        />
        <div className={styles.content}>
          <Scrollable>
            <div className={styles.padding}>
              <Divider/>
              <List dense className={styles.list}>
                {children}
              </List>
            </div>
          </Scrollable>
        </div>
      </div>
    </div>
    <div className={clsx(!collapsed && styles.hidden)}>
      <Actions className={styles.collapsedActions}>
        {collapsed && <Action
          tooltip={'Show filter menu'}
          onClick={() => {
            onCollapsedChange(false)
          }}
          className={styles.button}
        >
          <ArrowForwardIcon fontSize="small"/>
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

const useFilterSubMenusStyles = makeStyles(theme => {
  const widthSmall = 25
  const widthMedium = 32
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
    containerSmall: {
      '-webkit-transform': `translateX(${widthSmall}rem)`,
      transform: `translateX(${widthSmall}rem)`
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
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'column'
    },
    content: {
      flex: 1,
      minHeight: 0
    },
    menuSmall: {
      width: `${widthSmall}rem`
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
  const { useResults } = useSearchContext()
  const nResults = useResults()?.pagination?.total
  const styles = useFilterSubMenusStyles()
  const { open, onOpenChange, size, collapsed } = useContext(filterMenuContext)
  const [menuStyle, containerStyle] = {
    small: [styles.menuSmall, styles.containerSmall],
    medium: [styles.menuMedium, styles.containerMedium],
    large: [styles.menuLarge, styles.containerLarge]
  }[size]

  return <Paper
    elevation={4}
    className={clsx(styles.root, open && containerStyle)}
  >
    <div className={clsx(styles.menu, menuStyle, collapsed && styles.collapsed)}>
      <FilterMenuTopHeader
        title={isNil(nResults) ? 'loading...' : pluralize('result', nResults, true)}
        actions={<Action
          tooltip="Close submenu"
          onClick={() => { onOpenChange(false) }}
          className={styles.button}
        >
          <ArrowBackIcon fontSize="small"/>
        </Action>}
      />
      <div className={styles.content}>
        {children}
      </div>
    </div>
  </Paper>
})
FilterSubMenus.propTypes = {
  sizes: PropTypes.object,
  children: PropTypes.node
}

const useFilterSubMenuStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column'
  },
  hidden: {
    display: 'none'
  },
  padding: {
    paddingBottom: theme.spacing(1.5),
    paddingLeft: theme.spacing(paddingHorizontal),
    paddingRight: theme.spacing(paddingHorizontal)
  },
  content: {
    flex: 1,
    minHeight: 0
  }
}))

export const FilterSubMenu = React.memo(({
  value,
  size,
  actions,
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
    <FilterMenuHeader title={selected} actions={actions}/>
    <div className={styles.content}>
      <Scrollable>
        <div className={styles.padding}>
          {children}
        </div>
      </Scrollable>
    </div>
  </div>
})
FilterSubMenu.propTypes = {
  value: PropTypes.string,
  size: PropTypes.oneOf(['small', 'medium', 'large']),
  actions: PropTypes.node,
  children: PropTypes.node
}
FilterSubMenu.defaultProps = {
  size: 'small'
}

export default FilterSubMenu
