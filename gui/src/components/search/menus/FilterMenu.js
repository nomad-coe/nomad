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
import ReplayIcon from '@material-ui/icons/Replay'
import MoreVert from '@material-ui/icons/MoreVert'
import CodeIcon from '@material-ui/icons/Code'
import Scrollable from '../../visualization/Scrollable'
import FilterSummary from '../FilterSummary'
import FilterSettings from './FilterSettings'
import { Actions, ActionHeader, Action } from '../../Actions'
import { useSearchContext } from '../SearchContext'
import { pluralize } from '../../../utils'
import { isNil } from 'lodash'
import { SourceApiCall, SourceApiDialogButton, SourceDialogDivider, SourceJsonCode } from '../../buttons/SourceDialogButton'

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
      paddingTop: theme.spacing(1.2),
      paddingRight: theme.spacing(paddingHorizontal),
      paddingLeft: theme.spacing(paddingHorizontal)
    },
    title: {
      display: 'flex',
      alignItems: 'center',
      fontSize: '0.75rem',
      marginTop: -3,
      marginBottom: -3
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
      paddingBottom: theme.spacing(1.1),
      paddingTop: theme.spacing(0.3),
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
  className,
  'data-testid': testID
}) => {
  const styles = useFilterMenuHeaderStyles()
  return <div className={clsx(className, styles.root)} data-testid={testID}>
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
  className: PropTypes.string,
  'data-testid': PropTypes.string
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
  const [size, setSize] = useState('m')

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
  const { useResetFilters, useRefresh, useApiData } = useSearchContext()
  const styles = useFilterMenuItemsStyles()
  const { open, onOpenChange, collapsed, onCollapsedChange } = useContext(filterMenuContext)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const isSettingsOpen = Boolean(anchorEl)
  const resetFilters = useResetFilters()
  const refresh = useRefresh()
  const apiData = useApiData()

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
        <FilterMenuTopHeader/>
        <FilterMenuHeader
          title="Filters"
          actions={<>
            <Action
              tooltip="Refresh results"
              onClick={() => refresh()}
            >
              <ReplayIcon fontSize="small"/>
            </Action>
            <Action
              tooltip="Clear filters"
              onClick={() => resetFilters()}
            >
              <ClearIcon fontSize="small"/>
            </Action>
            <Action
              tooltip=""
              ButtonComponent={SourceApiDialogButton}
              ButtonProps={{
                tooltip: "API",
                maxWidth: "lg",
                fullWidth: true,
                icon: <CodeIcon fontSize="small"/>,
                ButtonProps: {
                  size: "small"
                }
              }}
            >
              <Typography>
                NOMAD uses the same query format throughout its API. This is the query
                based on the current filters:
              </Typography>
              <SourceJsonCode data={{owner: apiData?.body?.owner, query: apiData?.body?.query}}/>
              <SourceDialogDivider/>
              <Typography>
                One application of the above query is this API call. This is what is currently
                used to render this page and includes all displayed statistics data
                (aggregations).
              </Typography>
              <SourceApiCall
                {...apiData}
              />
            </Action>
            <Action
              tooltip={'Hide filter menu'}
              onClick={() => {
                onCollapsedChange(old => !old)
                onOpenChange(false)
              }}
              // className={styles.button}
            >
              <ArrowBackIcon fontSize="small"/>
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

const levelIndent = 1.8
const leftGutter = 3.0
const rightGutter = 2.35
const itemHeight = 2.6
const useFilterMenuItemStyles = makeStyles(theme => {
  return {
    root: {
      minHeight: `${itemHeight}rem`
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
      paddingLeft: theme.spacing(leftGutter),
      paddingRight: theme.spacing(rightGutter)
    },
    listItem: {
      position: 'relative',
      height: `${itemHeight}rem`
    },
    actions: {
      display: 'flex',
      flexDirection: 'column',
      paddingLeft: theme.spacing(leftGutter),
      paddingRight: theme.spacing(rightGutter)
    },
    divider: {
      width: '100%',
      backgroundColor: theme.palette.grey[300]
    }
  }
})
export const FilterMenuItem = React.memo(({
  id,
  label,
  group,
  onClick,
  actions,
  disableButton,
  level
}) => {
  const styles = useFilterMenuItemStyles()
  const theme = useTheme()
  const {filterGroups} = useSearchContext()
  const groupFinal = group || filterGroups[id]
  const { selected, open, onChange } = useContext(filterMenuContext)
  const handleClick = disableButton ? undefined : (onClick || onChange)
  const opened = open && id === selected

  return <div className={styles.root}>
    {(label || handleClick) &&
    <ListItem
      button={!!handleClick}
      className={styles.listItem}
      classes={{gutters: styles.gutters}}
      onClick={handleClick && (() => handleClick(id))}
    >
      {label && <ListItemText
        style={{marginLeft: theme.spacing(level * levelIndent)}}
        primaryTypographyProps={{
          color: opened ? 'primary' : 'initial',
          className: styles.label
        }}
        primary={label}
        data-testid={`menu-item-label-${id}`}
      />}
      {handleClick && <ListItemIcon className={styles.listIcon}>
        <NavigateNextIcon color={opened ? 'primary' : 'action'} className={styles.arrow}/>
      </ListItemIcon>}
    </ListItem>}
    {actions && <div
      className={styles.actions}
      style={{marginLeft: theme.spacing(level * levelIndent)}}
    >
      {actions}
    </div>}
    {groupFinal && <FilterSummary quantities={groupFinal}/>}
    <Divider className={styles.divider}/>
  </div>
})

FilterMenuItem.propTypes = {
  id: PropTypes.string,
  label: PropTypes.string,
  group: PropTypes.string,
  onClick: PropTypes.func,
  actions: PropTypes.node,
  disableButton: PropTypes.bool,
  level: PropTypes.number
}
FilterMenuItem.defaultProps = {
  level: 0
}

const useFilterSubMenusStyles = makeStyles(theme => {
  const widthS = 25
  const widthM = 32
  const widthL = 40
  const widthXL = 48
  return {
    root: {
      zIndex: 2,
      display: 'flex',
      flexDirection: 'column',
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      width: `${widthXL}rem`,
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
    containerS: {
      '-webkit-transform': `translateX(${widthS}rem)`,
      transform: `translateX(${widthS}rem)`
    },
    containerM: {
      '-webkit-transform': `translateX(${widthM}rem)`,
      transform: `translateX(${widthM}rem)`
    },
    containerL: {
      '-webkit-transform': `translateX(${widthL}rem)`,
      transform: `translateX(${widthL}rem)`
    },
    containerXL: {
      '-webkit-transform': `translateX(${widthXL}rem)`,
      transform: `translateX(${widthXL}rem)`
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
    menuS: {
      width: `${widthS}rem`
    },
    menuM: {
      width: `${widthM}rem`
    },
    menuL: {
      width: `${widthL}rem`
    },
    menuXL: {
      width: `${widthXL}rem`
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
    s: [styles.menuS, styles.containerS],
    m: [styles.menuM, styles.containerM],
    l: [styles.menuL, styles.containerL],
    xl: [styles.menuXL, styles.containerXL]
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
  id,
  label,
  size,
  actions,
  children
}) => {
  const styles = useFilterSubMenuStyles()
  const { selected, onSizeChange } = useContext(filterMenuContext)
  const visible = id === selected
  useEffect(() => {
    if (visible) {
      onSizeChange(size)
    }
  }, [size, visible, onSizeChange])

  return <div className={clsx(styles.root, !visible && styles.hidden)}>
    <FilterMenuHeader title={label} actions={actions} data-testid={`filter-menu-header-${id}`}/>
    <div className={styles.content}>
      <Scrollable>
        {children}
      </Scrollable>
    </div>
  </div>
})
FilterSubMenu.propTypes = {
  id: PropTypes.string,
  label: PropTypes.string,
  size: PropTypes.oneOf(['s', 'm', 'l', 'xl']),
  actions: PropTypes.node,
  children: PropTypes.node
}
FilterSubMenu.defaultProps = {
  size: 's'
}

export default FilterSubMenu
