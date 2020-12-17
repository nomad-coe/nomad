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

import React, { useContext, useState } from 'react'
import PropTypes from 'prop-types'
import { Button, makeStyles, MenuItem as MuiMenuItem, Menu as MuiMenu, ClickAwayListener, ListItemText } from '@material-ui/core'
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown'
import { useHistory } from 'react-router-dom'

const capitalize = (s) => {
  if (typeof s !== 'string') return ''
  return s.charAt(0).toUpperCase() + s.slice(1)
}

const MenuBarContext = React.createContext({})

const useMenuBarStyles = makeStyles(theme => ({
  root: {
    margin: theme.spacing(1),
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    flexWrap: 'nowrap',
    justifyContent: 'left'
  }
}))

export function MenuBar({children, selected}) {
  const classes = useMenuBarStyles()

  const [openMenu, setOpenMenu] = useState(null)
  const [openMenuAnchorEl, setOpenMenuAnchorEl] = useState(null)

  const handleClickMenu = (menuName, event) => {
    setOpenMenu(menuName)
    setOpenMenuAnchorEl(event.currentTarget)
  }

  const handleCloseMenus = () => {
    setOpenMenu(null)
    setOpenMenuAnchorEl(null)
  }

  const [selectedMenu, selectedMenuItem] = (selected || '').split('/')

  const context = {
    selectedMenu: selectedMenu,
    selectedMenuItem: selectedMenuItem,
    openMenu: openMenu,
    openMenuAnchorEl: openMenuAnchorEl,
    onClickMenu: handleClickMenu,
    onCloseMenu: handleCloseMenus
  }
  return <MenuBarContext.Provider value={context}>
    <div className={classes.root}>
      {children}
    </div>
  </MenuBarContext.Provider>
}
MenuBar.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  selected: PropTypes.string.isRequired
}

const useMenuBarItemStyles = makeStyles(theme => ({
  selected: {
    color: theme.palette.primary.main
  },
  link: {
    color: 'inherit',
    textDecoration: 'none'
  }
}))

export const MenuBarItem = React.forwardRef(({name, label, tooltip, route, href, onClick}, ref) => {
  const classes = useMenuBarItemStyles()

  label = label || capitalize(name)
  const context = useContext(MenuBarContext)
  const selected = context.selectedMenuItem === name

  const history = useHistory()

  const handleClick = () => {
    context.onCloseMenu(context.selectedMenu)
    if (route) {
      history.push(route)
    } else if (!href) {
      onClick()
    }
  }

  const item = <MuiMenuItem
    ref={ref}
    dense
    classes={{root: selected ? classes.selected : undefined}}
    onClick={handleClick}
  >
    <ListItemText primary={label} secondary={tooltip} />
  </MuiMenuItem>

  return href ? <a href={href} className={classes.link}>{item}</a> : item
})
MenuBarItem.propTypes = {
  name: PropTypes.string.isRequired,
  label: PropTypes.string,
  tooltip: PropTypes.string,
  icon: PropTypes.node,
  route: PropTypes.string,
  href: PropTypes.string,
  onClick: PropTypes.func
}

const useMenuBarMenuStyles = makeStyles(theme => ({
  root: {
    marginLeft: theme.spacing(1),
    '&:first-child': {
      marginLeft: theme.spacing(0)
    }
  },
  menuPopover: {
    pointerEvents: 'none'
  },
  menuPaper: {
    pointerEvents: 'all'
  }
}))

export function MenuBarMenu({name, label, children}) {
  const classes = useMenuBarMenuStyles()
  label = label || capitalize(name)
  const context = useContext(MenuBarContext)

  const selected = context.selectedMenu === name
  const open = context.openMenu === name
  const aMenuIsOpen = context.openMenu !== null
  const anchorEl = open ? context.openMenuAnchorEl : undefined

  const handleClick = (event) => {
    context.onClickMenu(name, event)
  }

  const handleMouseOver = (event) => {
    if (aMenuIsOpen) {
      context.onClickMenu(name, event)
    }
  }

  const handleClose = () => {
    if (open) {
      context.onCloseMenu(name)
    }
  }

  return <React.Fragment>
    <ClickAwayListener onClickAway={handleClose}>
      <Button
        className={classes.root}
        onClick={handleClick}
        onMouseOver={handleMouseOver}
        disableElevation
        size="small"
        color={selected ? 'primary' : undefined}
        endIcon={<KeyboardArrowDownIcon />}
      >
        {label}
      </Button>
    </ClickAwayListener>
    <MuiMenu
      PopoverClasses={{root: classes.menuPopover, paper: classes.menuPaper}}
      elevation={1}
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{vertical: 'bottom', horizontal: 'left'}}
      keepMounted
      open={open}
      onClose={handleClose}
    >
      {children}
    </MuiMenu>
  </React.Fragment>
}
MenuBarMenu.propTypes = {
  name: PropTypes.string.isRequired,
  label: PropTypes.string,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}
