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
import PropTypes from 'prop-types'
import {
  Button,
  Box,
  makeStyles,
  MenuList,
  MenuItem,
  Menu,
  ListItemText,
  Typography
} from '@material-ui/core'
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown'
import { matchPath, useHistory, useLocation } from 'react-router-dom'

const useMenuBarStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    flexWrap: 'nowrap',
    justifyContent: 'left',
    boxSizing: 'border-box'
  }
}))

/**
 * Reusable menubar e.g. for navigation.
*/
export function MenuBar({children}) {
  const classes = useMenuBarStyles()

  return <div className={classes.root}>
    {children}
  </div>
}
MenuBar.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export const MenuBarList = React.memo(({header, children}) => {
  return <MenuList
      subheader={
        <Box marginLeft={1.5}>
          <Typography variant="button">{header}</Typography>
        </Box>
      }
    >
      {children}
    </MenuList>
})
MenuBarList.propTypes = {
  header: PropTypes.string,
  children: PropTypes.node
}

const useMenuBarItemStyles = makeStyles(theme => ({
  selected: {
    color: theme.palette.primary.main
  }
}))

export const MenuBarItem = React.forwardRef(({label, tooltip, route, href}, ref) => {
  const classes = useMenuBarItemStyles()
  const {pathname} = useLocation()
  const selected = matchPath(pathname, route) && true

  const history = useHistory()

  const handleClick = () => {
    if (route) {
      history.push(route)
    }
  }

  return <MenuItem
    data-testid={label}
    ref={ref}
    component={href ? 'a' : 'li'}
    href={href}
    dense
    classes={{root: selected ? classes.selected : undefined}}
    onClick={handleClick}
  >
    <ListItemText primary={label} secondary={tooltip} />
  </MenuItem>
})
MenuBarItem.propTypes = {
  label: PropTypes.string.isRequired,
  tooltip: PropTypes.string,
  icon: PropTypes.node,
  route: PropTypes.string,
  href: PropTypes.string
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

export function MenuBarMenu({label, children, route}) {
  const classes = useMenuBarMenuStyles()
  const {pathname} = useLocation()
  const [anchorEl, setAnchorEl] = useState(null)
  const selected = matchPath(pathname, route) && true

  const handleOpen = (event) => {
    setAnchorEl(event.currentTarget)
  }

  const handleClose = () => {
    setAnchorEl(null)
  }

  return <div onMouseLeave={handleClose}>
    <Button
      className={classes.root}
      onMouseEnter={handleOpen}
      disableElevation
      size="small"
      color={selected ? 'primary' : undefined}
      endIcon={<KeyboardArrowDownIcon />}
    >
      {label}
    </Button>
    <Menu
      data-testid={label}
      PopoverClasses={{root: classes.menuPopover, paper: classes.menuPaper}}
      elevation={2}
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{vertical: 'bottom', horizontal: 'left'}}
      keepMounted
      open={Boolean(anchorEl)}
      onClose={handleClose}
      onClick={handleClose}
    >
      {children}
    </Menu>
  </div>
}
MenuBarMenu.propTypes = {
  label: PropTypes.string.isRequired,
  route: PropTypes.string.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}
