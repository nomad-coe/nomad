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
import React, { useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { Alert } from '@material-ui/lab'
import { has } from 'lodash'
import {
  FilterMenu,
  FilterMenuItem,
  FilterMenuItems,
  FilterSubMenus
} from './FilterMenu'
import { makeStyles } from '@material-ui/core/styles'
import FilterSubMenuElements from './FilterSubMenuElements'
import FilterSubMenuStructure from './FilterSubMenuStructure'
import FilterSubMenuMethod from './FilterSubMenuMethod'
import FilterSubMenuPrecision from './FilterSubMenuPrecision'
import FilterSubMenuDFT from './FilterSubMenuDFT'
import FilterSubMenuTB from './FilterSubMenuTB'
import FilterSubMenuGW from './FilterSubMenuGW'
import FilterSubMenuBSE from './FilterSubMenuBSE'
import FilterSubMenuDMFT from './FilterSubMenuDMFT'
import FilterSubMenuEELS from './FilterSubMenuEELS'
import FilterSubMenuElectronic from './FilterSubMenuElectronic'
import FilterSubMenuSolarCell from './FilterSubMenuSolarCell'
import FilterSubMenuCatalyst from './FilterSubMenuCatalystProperties'
import FilterSubMenuVibrational from './FilterSubMenuVibrational'
import FilterSubMenuMechanical from './FilterSubMenuMechanical'
import FilterSubMenuMolecularDynamics from './FilterSubMenuMolecularDynamics'
import FilterSubMenuELN from './FilterSubMenuELN'
import FilterSubMenuAuthor from './FilterSubMenuAuthor'
import FilterSubMenuMetadata from './FilterSubMenuMetadata'
import FilterSubMenuOptimade from './FilterSubMenuOptimade'
import { useSearchContext } from '../SearchContext'
import { delay } from '../../../utils'
import FilterSubMenuGeometryOptimization from './FilterSubMenuGeometryOptimization'
import InputCheckbox from '../input/InputCheckbox'
import FilterSubMenuCustomQuantities from './FilterSubMenuCustomQuantities'

export const menuMap = {
  elements: FilterSubMenuElements,
  structure: FilterSubMenuStructure,
  method: FilterSubMenuMethod,
  precision: FilterSubMenuPrecision,
  dft: FilterSubMenuDFT,
  tb: FilterSubMenuTB,
  gw: FilterSubMenuGW,
  bse: FilterSubMenuBSE,
  dmft: FilterSubMenuDMFT,
  eels: FilterSubMenuEELS,
  electronic: FilterSubMenuElectronic,
  solarcell: FilterSubMenuSolarCell,
  heterogeneouscatalyst: FilterSubMenuCatalyst,
  vibrational: FilterSubMenuVibrational,
  mechanical: FilterSubMenuMechanical,
  molecular_dynamics: FilterSubMenuMolecularDynamics,
  geometry_optimization: FilterSubMenuGeometryOptimization,
  eln: FilterSubMenuELN,
  custom_quantities: FilterSubMenuCustomQuantities,
  author: FilterSubMenuAuthor,
  metadata: FilterSubMenuMetadata,
  optimade: FilterSubMenuOptimade
}

const useFilterMainMenuStyles = makeStyles(theme => ({
  combine: {
  }
}))

/**
 * Swipable menu that shows the available filters on the left side of the
 * screen.
 */
const FilterMainMenu = React.memo(({
  open,
  onOpenChange,
  collapsed,
  onCollapsedChange
}) => {
  const [value, setValue] = React.useState()
  const {filterMenus} = useSearchContext()
  const [loaded, setLoaded] = useState(false)
  const styles = useFilterMainMenuStyles()

  // Rendering the submenus is delayed on the event queue: this makes loading
  // the search page more responsive by first loading everything else.
  useEffect(() => {
    delay(() => { setLoaded(true) })
  }, [])

  // The shown menu items
  const menuItems = useMemo(() => {
    return filterMenus?.options
     ? Object.values(filterMenus.options).map(option => {
        return <FilterMenuItem
          key={option.key}
          id={option.key}
          label={option.label}
          level={option.level}
          disableButton={!has(menuMap, option.key)}
          actions={option?.actions?.options && Object.values(option.actions.options)
            .map((action) => {
              const content = action.type === 'checkbox'
                ? <InputCheckbox
                  key={action.key}
                  quantity={action.quantity}
                  description={action.tooltip}
                  label={action.label}
                  className={styles.combine}
                ></InputCheckbox>
                : null
              return content
          })}
        />
      })
    : <Alert severity="warning">
      No search menus defined within this search context. Ensure that all GUI artifacts are created.
    </Alert>
  }, [filterMenus, styles])

  // The shown submenus
  const subMenus = useMemo(() => {
    return filterMenus?.options
      ? Object.values(filterMenus.options)
        .filter(option => menuMap[option.key])
        .map(option => {
          const Comp = menuMap[option.key]
          return <Comp
            key={option.key}
            id={option.key}
            label={option.label}
            size={option.size}
          />
        })
      : null
  }, [filterMenus])

  return <FilterMenu
    selected={value}
    onSelectedChange={setValue}
    open={open}
    onOpenChange={onOpenChange}
    collapsed={collapsed}
    onCollapsedChange={onCollapsedChange}
  >
    <FilterMenuItems>
      {menuItems}
    </FilterMenuItems>
    <FilterSubMenus>
      {loaded && subMenus}
    </FilterSubMenus>
  </FilterMenu>
})
FilterMainMenu.propTypes = {
  open: PropTypes.bool,
  onOpenChange: PropTypes.func,
  collapsed: PropTypes.bool,
  onCollapsedChange: PropTypes.func
}

export default FilterMainMenu
