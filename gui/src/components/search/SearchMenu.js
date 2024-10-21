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
import React, { useEffect, useMemo, useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import { isEmpty } from 'lodash'
import {
  Menu as MenuMUI
} from '@material-ui/core'
import {
  Menu,
  MenuHeader,
  MenuItem,
  MenuContent,
  MenuSubMenus,
  MenuSettings
} from './menus/Menu'
import ArrowBackIcon from '@material-ui/icons/ArrowBack'
import MoreVert from '@material-ui/icons/MoreVert'
import { Action } from '../Actions'
import { useSearchContext } from './SearchContext'
import { delay, camelCase } from '../../utils'
import InputTerms from './input/InputTerms'
import InputHistogram from './input/InputHistogram'
import InputPeriodicTable from './input/InputPeriodicTable'
import InputNestedObject from './input/InputNestedObject'
import InputVisibility from './input/InputVisibility'
import InputDefinitions from './input/InputDefinitions'
import InputOptimade from './input/InputOptimade'
import InputCustomQuantities from './input/InputCustomQuantities'
import { InputGrid, InputGridItem } from './input/InputGrid'

function createItems(items, visible) {
  return (items || []).map((item, index) => {
    let props = Object.fromEntries(Object.entries(item).map(([key, value]) => [camelCase(key), value]))
    props = {visible, ...props}
    let Comp
    switch (item.type) {
      case 'terms':
        Comp = InputTerms
        break
      case 'histogram':
        Comp = InputHistogram
        break
      case 'periodic_table':
        Comp = InputPeriodicTable
        break
      case 'nested_object':
        Comp = InputNestedObject
        props.children = createItems(item.items, visible)
        break
      case 'visibility':
        Comp = InputVisibility
        break
      case 'definitions':
        Comp = InputDefinitions
        break
      case 'optimade':
        Comp = InputOptimade
        break
      case 'custom_quantities':
        Comp = InputCustomQuantities
        break
      case 'menu':
        Comp = MenuItem
        props.disableButton = isEmpty(item.items)
        props.id = index
        props.level = item.indentation
        break
      default:
        throw Error(`Unknown menu item type: ${item.type}.`)
    }

    return <InputGridItem
      key={index}
      disablePadding={item.type === 'menu'}
      xs={item.type === 'menu' ? 12 : item.width}>
        {Comp && <Comp {...props}/>}
    </InputGridItem>
  })
}
/**
 * Creates an instance of FilterMenu based on the menu in the SearchContext.
 */
const SearchMenu = React.memo(() => {
  const [selected, setSelected] = React.useState()
  const [isSubMenuOpen, setIsSubMenuOpen] = React.useState(false)
  const [loaded, setLoaded] = useState(false)
  const {menu, useSetIsMenuOpen, useIsMenuOpen, useIsCollapsed, useSetIsCollapsed} = useSearchContext()
  const [isMenuOpen, setIsMenuOpen] = [useIsMenuOpen(), useSetIsMenuOpen()]
  const [isCollapsed, setIsCollapsed] = [useIsCollapsed(), useSetIsCollapsed()]
  const handleMenuCollapse = useCallback(() => setIsCollapsed(old => !old), [setIsCollapsed])

  // Rendering the submenus is delayed on the event queue: this makes loading
  // the search page more responsive by first loading everything else.
  useEffect(() => {
    delay(() => { setLoaded(true) })
  }, [])

  useEffect(() => {
    setIsSubMenuOpen(isMenuOpen)
  }, [isMenuOpen])

  const [anchorEl, setAnchorEl] = React.useState(null)
  const isSettingsOpen = Boolean(anchorEl)

  // Callbacks
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])
  const handleOpenChange = useCallback((value) => {
    setIsMenuOpen(value)
    setIsSubMenuOpen(value)
  }, [setIsMenuOpen])

  // Create the list of menu items
  const menuItems = useMemo(
    () => createItems(menu?.items || [], true),
    [menu]
  )

  // Create the list of submenus
  const subMenus = useMemo(() => {
    return (menu?.items || [])
      .map((item, index) =>
        item.type === 'menu'
          ? <SearchSubMenu
            key={index}
            menu={item}
            selected={selected}
            open={isSubMenuOpen}
            visible={index === selected}
            onOpenChange={handleOpenChange}
          />
          : null
      )
  }, [menu, selected, handleOpenChange, isSubMenuOpen])

  return <Menu
    size={menu?.size}
    open={true}
    collapsed={isCollapsed}
    onCollapsedChanged={setIsCollapsed}
    subMenuOpen={isSubMenuOpen}
    visible={true}
    selected={selected}
    onSelectedChange={(value) => {
      setIsMenuOpen(old => value !== selected ? true : !old)
      setSelected(value)
      setIsSubMenuOpen(old => value !== selected ? true : !old)
    }}
  >
    <MenuHeader
      title="Filters"
      actions={<>
        <Action
          tooltip={'Collapse menu'}
          onClick={handleMenuCollapse}
        >
          <ArrowBackIcon fontSize="small"/>
        </Action>
        <Action
          tooltip="Options"
          onClick={openMenu}
        >
          <MoreVert fontSize="small"/>
        </Action>
        <MenuMUI
          anchorEl={anchorEl}
          open={isSettingsOpen}
          onClose={closeMenu}
          getContentAnchorEl={null}
          anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
          transformOrigin={{ vertical: 'top', horizontal: 'right' }}
          keepMounted
        >
          <div>
            <MenuSettings/>
          </div>
        </MenuMUI>
      </>}
    />
    <MenuContent>
      <InputGrid>{menuItems}</InputGrid>
    </MenuContent>
    <MenuSubMenus>
      {loaded && subMenus}
    </MenuSubMenus>
  </Menu>
})

/**
 * Submenu that pops on the right side of the parent menu.
 */
const SearchSubMenu = React.memo(({
  menu,
  open,
  visible,
  onOpenChange
}) => {
  const [selected, setSelected] = React.useState()
  const [isSubMenuOpen, setIsSubMenuOpen] = React.useState(false)

  // Callbacks
  const handleClose = useCallback(() => onOpenChange(false), [onOpenChange])
  const handleOpenChange = useCallback((value) => {
    setIsSubMenuOpen(value)
  }, [])

  // Create the list of menu items
  const menuItems = useMemo(
    () => createItems(menu.items, visible),
    [menu.items, visible]
  )

  // Create the list of submenus
  const subMenus = useMemo(() => {
    return (menu.items || [])
      .map((item, index) =>
        item.type === 'menu'
          ? <SearchSubMenu
            key={index}
            menu={item}
            open={open && isSubMenuOpen}
            visible={visible && index === selected}
            onOpenChange={handleOpenChange}
          />
          : null
      )
  }, [menu, selected, handleOpenChange, isSubMenuOpen, open, visible])

  return <Menu
    size={menu.size}
    open={open}
    subMenuOpen={open && isSubMenuOpen}
    visible={visible}
    selected={selected}
    onSelectedChange={(value) => {
      setSelected(value)
      setIsSubMenuOpen(old => value !== selected ? true : !old)
    }}
  >
    <MenuHeader
      title={menu.title}
      actions={<>
        <Action
          tooltip={'Close menu'}
          onClick={handleClose}
        >
          <ArrowBackIcon fontSize="small"/>
        </Action>
      </>}
    />
    <MenuContent>
      <InputGrid>{menuItems}</InputGrid>
    </MenuContent>
    <MenuSubMenus>
      {subMenus}
    </MenuSubMenus>
  </Menu>
})
SearchSubMenu.propTypes = {
  menu: PropTypes.object,
  open: PropTypes.bool,
  visible: PropTypes.bool,
  selected: PropTypes.number,
  onOpenChange: PropTypes.func
}

export default SearchMenu
