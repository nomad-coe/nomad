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
import React, { useState, useEffect, isValidElement, Children, useCallback, useContext } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Paper,
  Typography,
  Checkbox,
  MenuItem as MenuItemMUI,
  FormControlLabel
} from '@material-ui/core'
import { delay } from '../../../utils'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import ArrowForwardIcon from '@material-ui/icons/ArrowForward'
import Scrollable from '../../visualization/Scrollable'
import { Action, Actions, ActionHeader } from '../../Actions'
import FilterTitle from '../FilterTitle'
import { useSearchContext } from '../SearchContext'

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

export const menuContext = React.createContext()
export const paddingHorizontal = 1.5
const collapseWidth = '3rem'

/**
 * Provides a context for a menu.
 */
const useMenuStyles = (props) => makeStyles(theme => {
  return {
    root: {
      boxSizing: 'border-box',
      display: 'flex',
      height: '100%',
      flexDirection: 'column',
      '-webkit-transform': `translateX(-${props.width})`,
      'transform': `translateX(-${props.width})`,
      width: `${props.width}`,
      transition: `transform 225ms`,
      willChange: 'transform',
      visibility: 'hidden'
    },
    menu: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'column',
      marginTop: '1px',
      width: `${props.width}`
    },
    open: {
      '-webkit-transform': 'none',
      transform: 'none'
    },
    container: {
      height: '100%'
    },
    containerCollapsed: {
      width: collapseWidth
    },
    collapsed: {
      '-webkit-transform': `translateX(calc(-${props.width} + ${collapseWidth}))`,
      'transform': `translateX(calc(-${props.width} + 50px))`
    },
    visible: {
      visibility: 'visible'
    }
  }
})
export const Menu = React.memo(({
  size,
  open,
  collapsed,
  onCollapsedChanged,
  subMenuOpen,
  visible,
  onOpenChange,
  selected,
  onSelectedChange,
  className,
  children
}) => {
  const width = {
    'xs': '17rem',
    'sm': '21rem',
    'md': '25rem',
    'lg': '29rem',
    'xl': '33rem',
    'xxl': '45rem'
  }[size] || size || '21rem'
  const styles = useMenuStyles({width})()

  return <div className={clsx(styles.container, collapsed && styles.containerCollapsed)}>
    <Paper
      elevation={open ? 4 : 0}
      className={clsx(className, styles.root, open && styles.open, collapsed && styles.collapsed, visible && styles.visible)}
    >
      <menuContext.Provider value={{
        collapsed,
        selected,
        setSelected: onSelectedChange,
        setCollapsed: onCollapsedChanged,
        subMenuOpen,
        open,
        setOpen: onOpenChange,
        size
      }}>
      {children}
      </menuContext.Provider>
    </Paper>
  </div>
})
Menu.propTypes = {
  size: PropTypes.string,
  open: PropTypes.bool,
  collapsed: PropTypes.bool,
  onCollapsedChanged: PropTypes.func,
  subMenuOpen: PropTypes.bool,
  visible: PropTypes.bool,
  onOpenChange: PropTypes.func,
  selected: PropTypes.number,
  onSelectedChange: PropTypes.func,
  className: PropTypes.string,
  children: PropTypes.node
}

/**
 * Header for filter menus. Contains panel actions and an overline text.
 */
const useMenuHeaderStyles = makeStyles(theme => {
  return {
    root: {
      height: theme.spacing(3),
      paddingBottom: theme.spacing(1.2),
      paddingTop: theme.spacing(1.4),
      paddingLeft: theme.spacing(paddingHorizontal),
      paddingRight: theme.spacing(paddingHorizontal),
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'center',
      zIndex: 4,
      backgroundColor: theme.palette.background.paper,
      boxShadow: `1px 0px 0px 0px #e0e0e0`
    },
    title: {
      display: 'flex',
      alignItems: 'center',
      fontSize: '0.90rem'
    }
  }
})
export const MenuHeader = React.memo(({
  title,
  actions,
  className,
  'data-testid': testID
}) => {
  const styles = useMenuHeaderStyles()
  const { collapsed, setCollapsed } = useContext(menuContext)

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <Actions>
      <ActionHeader>
        <Typography className={styles.title} variant="button">{title}</Typography>
      </ActionHeader>
      {collapsed
        ? <Action
          tooltip={'Show menu'}
          onClick={() => setCollapsed(false)}
        >
          <ArrowForwardIcon fontSize="small"/>
        </Action>
        : actions
      }
    </Actions>
  </div>
})

MenuHeader.propTypes = {
  overlineTitle: PropTypes.string,
  title: PropTypes.string,
  topAction: PropTypes.node,
  actions: PropTypes.node,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

/**
 * Menu content that is wrapped in a customized scrollable area.
 */
const useMenuContentStyles = makeStyles(theme => {
  return {
    root: {
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'column',
      width: '100%',
      height: '100%',
      position: 'relative',
      backgroundColor: theme.palette.background.paper,
      zIndex: 3
    },
    headerTextVertical: {
      display: 'flex',
      alignItems: 'center',
      position: 'absolute',
      right: '0.07rem',
      top: '2.5rem',
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
      boxShadow: `1px 0px 0px 0px #e0e0e0`
    },
    button: {
      marginRight: 0,
      '-webkit-transform': 'none',
      transform: 'none',
      transition: 'transform 250ms',
      willChange: 'transform'
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
    hidden: {
      visibility: 'hidden'
    },
    collapsedText: {
      display: 'flex',
      alignItems: 'center',
      position: 'absolute',
      right: '0.07rem',
      top: 0,
      height: collapseWidth,
      transform: 'rotate(90deg)'
    }

  }
})
export const MenuContent = React.memo(({
  className,
  children
}) => {
  const styles = useMenuContentStyles()
  const { open, collapsed } = useContext(menuContext)
  // Unfortunately the ClickAwayListener does not play nicely together with
  // Menus/Select/Popper. When using Portals, the clicks are registered wrong.
  // When Portals are disabled (disablePortal), their positioning goes haywire.
  // The clicks outside are thus detected by individual event listeners that
  // toggle the menu state.
  return <div className={clsx(className, styles.root)}>
    <div className={clsx(styles.menu, open && styles.menuBorder)}>
      <div className={styles.container}>
        <div className={clsx(styles.content, collapsed && styles.hidden)}>
          <Scrollable>
            <List dense className={styles.list}>
              {children}
            </List>
          </Scrollable>
        </div>
      {collapsed && <Typography
          className={clsx(styles.collapsedText)}
          variant="button"
        >Filters
        </Typography>
      }
      </div>
    </div>
  </div>
})
MenuContent.propTypes = {
  className: PropTypes.string,
  children: PropTypes.node
}

/**
 * Menu submenus. Ensures that submenus are correctly placed underneath the
 * menu and their loading is delayed.
 */
const useMenuSubMenusStyles = makeStyles(theme => {
  return {
    root: {
      position: 'absolute',
      top: '0',
      bottom: '0',
      left: '100%'
    },
    relative: {
      position: 'relative',
      width: '100%',
      height: '100%'
    },
    absolute: {
      position: 'absolute',
      top: 0,
      bottom: 0,
      right: 0,
      left: 0
    }
  }
})
export const MenuSubMenus = React.memo(({
  className,
  children
}) => {
  const styles = useMenuSubMenusStyles()
  const [loaded, setLoaded] = useState(false)

  // Rendering the submenus is delayed on the event queue: this makes loading
  // the search page more responsive by first loading everything else.
  useEffect(() => {
    delay(() => { setLoaded(true) })
  }, [])

  return loaded
    ? <div className={clsx(className, styles.root)}>
    <div className={styles.relatives}>
      {Children.map(children, (child) => {
          if (isValidElement(child)) {
            return (
              <div className={styles.absolute}>
                <child.type {...child.props}/>
              </div>
            )
          }
      })}
    </div>
  </div>
  : null
})
MenuSubMenus.propTypes = {
  className: PropTypes.string,
  children: PropTypes.node
}

/**
 * Menu item.
 */
const levelIndent = 1.8
const useMenuItemStyles = makeStyles(theme => {
  return {
    root: {
    },
    listIcon: {
      fontsize: '1rem',
      minWidth: '1.5rem'
    },
    arrow: {
      marginLeft: theme.spacing(1),
      fontSize: '1.5rem'
    },
    listItem: {
      paddingTop: theme.spacing(0.75),
      paddingBottom: theme.spacing(0.75),
      position: 'relative'
    },
    actions: {
      display: 'flex',
      flexDirection: 'column'
    }
  }
})
export const MenuItem = React.memo(({
  id,
  title,
  disableButton,
  level
}) => {
  const styles = useMenuItemStyles()
  const theme = useTheme()
  const { selected, setSelected, subMenuOpen } = useContext(menuContext)
  const opened = subMenuOpen && id === selected
  const handleClick = useCallback(() => {
    if (disableButton) return
    setSelected?.(id)
  }, [disableButton, id, setSelected])

  return <div className={styles.root}>
    {(title || !disableButton) &&
    <ListItem
      button={!disableButton}
      className={styles.listItem}
      onClick={handleClick}
    >
      {title && <ListItemText
        style={{marginLeft: theme.spacing(level * levelIndent)}}
        primaryTypographyProps={{
          color: opened ? 'primary' : 'initial'
        }}
        primary={<FilterTitle label={title}/>}
        data-testid={`menu-item-label-${id}`}
      />}
      {!disableButton && <ListItemIcon className={styles.listIcon}>
        <NavigateNextIcon color={opened ? 'primary' : 'action'} className={styles.arrow}/>
      </ListItemIcon>}
    </ListItem>}
  </div>
})

MenuItem.propTypes = {
  id: PropTypes.number,
  title: PropTypes.string,
  disableButton: PropTypes.bool,
  level: PropTypes.number
}
MenuItem.defaultProps = {
  level: 0
}

/**
 * Settings for a menu.
 */
export const MenuSettings = React.memo(() => {
  const {
    useIsStatisticsEnabled,
    useSetIsStatisticsEnabled
  } = useSearchContext()
  const [isStatisticsEnabled, setIsStatisticsEnabled] = [
    useIsStatisticsEnabled(),
    useSetIsStatisticsEnabled()
  ]

  const handleStatsChange = useCallback((event, value) => {
    setIsStatisticsEnabled(value)
  }, [setIsStatisticsEnabled])

  return <MenuItemMUI>
    <FormControlLabel
      control={<Checkbox
        checked={isStatisticsEnabled}
        onChange={handleStatsChange}
      />}
      label="Show advanced statistics"
    />
  </MenuItemMUI>
})
