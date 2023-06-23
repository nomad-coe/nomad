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

import { isEmpty } from 'lodash'
import React from 'react'
import PropTypes from 'prop-types'
import { Divider } from '@material-ui/core'
import { MenuBarItem, MenuBarList, MenuBarMenu } from './MenuBar'

const MenuBarRoute = React.memo(({
  menu,
  label,
  initialFilters,
  children
}) => {
    let filters = initialFilters && Object.entries(initialFilters).length > 0 && Object.entries(initialFilters).map(x => x.join('=')).join('&')
    filters = filters ? `?${filters}` : ''
    if (!menu.menu) return null
    // Gather all categories for the items in this menu
    const categories = {}
    menu.routes
      .filter(item => item.category)
      .forEach(item => {
        if (!categories[item.category]) categories[item.category] = []
        categories[item.category].push(item)
      })

    // Items without a category are laid out linearly
    let content
    if (isEmpty(categories)) {
      content = menu.routes.filter(route => route.menu).map((itemRoute, i) => (
        <MenuBarItem
          key={i} label={itemRoute.menu} tooltip={itemRoute.tooltip}
          route={itemRoute.path && `/${menu.path}/${itemRoute.path}${filters}`}
          href={itemRoute.href}
        />
      ))
    // Categorized items are laid out in columns
    } else {
      content = Object.entries(categories).map(([category, values], index) => {
          return <div key={category}>
            {index ? <Divider style={{marginBottom: 8}}/> : null}
            <MenuBarList header={category}>
              {values.map(item => (
                <MenuBarItem
                  key={item.path} label={item.menu} tooltip={item.tooltip}
                  route={item.path && `/${menu.path}/${item.path}${filters}`}
                  href={item.href}
                />
              ))}
            </MenuBarList>
          </div>
        })
    }

  return <MenuBarMenu key={menu.path} label={label} route={'/' + menu.path} icon={children}>
        {content}
      </MenuBarMenu>
})

MenuBarRoute.propTypes = {
  menu: PropTypes.object.isRequired,
  label: PropTypes.string.isRequired,
  // The 'initialFilters' property is optional and should be a object representing
  // the initial set of filters for the menu. For ex, it could be set to '{upload_id:abc}'
  // to preselect items based on the 'upload_id' parameter with the value 'abc'
  initialFilters: PropTypes.object,
  children: PropTypes.node
}

export default MenuBarRoute
