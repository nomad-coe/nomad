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
import React, { useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import {
  FilterMenu,
  FilterMenuItem,
  FilterMenuItems,
  FilterSubMenus
} from './FilterMenu'
import { makeStyles } from '@material-ui/core/styles'
import FilterSubMenuMaterial from './FilterSubMenuMaterial'
import FilterSubMenuElements from './FilterSubMenuElements'
import FilterSubMenuSymmetry from './FilterSubMenuSymmetry'
import FilterSubMenuMethod from './FilterSubMenuMethod'
import FilterSubMenuSimulation from './FilterSubMenuSimulation'
import FilterSubMenuDFT from './FilterSubMenuDFT'
import FilterSubMenuGW from './FilterSubMenuGW'
import FilterSubMenuEELS from './FilterSubMenuEELS'
import FilterSubMenuElectronic from './FilterSubMenuElectronic'
import FilterSubMenuVibrational from './FilterSubMenuVibrational'
import FilterSubMenuMechanical from './FilterSubMenuMechanical'
import FilterSubMenuSpectroscopy from './FilterSubMenuSpectroscopy.js'
import FilterSubMenuELN from './FilterSubMenuELN.js'
import FilterSubMenuAuthor from './FilterSubMenuAuthor'
import FilterSubMenuAccess from './FilterSubMenuAccess'
import FilterSubMenuDataset from './FilterSubMenuDataset'
import FilterSubMenuIDs from './FilterSubMenuIDs'
import FilterSubMenuArchive from './FilterSubMenuArchive'
import {
  labelMaterial,
  labelElements,
  labelSymmetry,
  labelMethod,
  labelSimulation,
  labelDFT,
  labelGW,
  labelExperiment,
  labelEELS,
  labelProperties,
  labelElectronic,
  labelVibrational,
  labelMechanical,
  labelELN,
  labelAuthor,
  labelDataset,
  labelIDs,
  labelAccess,
  labelSpectroscopy,
  labelArchive,
  labelGeometryOptimization
} from '../FilterRegistry'
import { useSearchContext } from '../SearchContext'
import InputCheckbox from '../input/InputCheckbox'
import { delay } from '../../../utils'
import FilterSubMenuGeometryOptimization from './FilterSubMenuGeometryOptimization'

/**
 * Swipable menu that shows the available filters on the left side of the
 * screen.
 */
const useStyles = makeStyles(theme => ({
  combine: {
    paddingLeft: theme.spacing(2)
  }
}))
const FilterMainMenu = React.memo(({
  open,
  onOpenChange,
  collapsed,
  onCollapsedChange
}) => {
  const [value, setValue] = React.useState()
  const {resource} = useSearchContext()
  const styles = useStyles()
  const [loaded, setLoaded] = useState(false)

  // Rendering the submenus is delayed: this makes loading the search page more
  // responsive by first loading everything else.
  useEffect(() => {
    delay(() => { setLoaded(true) })
  }, [])

  return <FilterMenu
    selected={value}
    onSelectedChange={setValue}
    open={open}
    onOpenChange={onOpenChange}
    collapsed={collapsed}
    onCollapsedChange={onCollapsedChange}
  >
    <FilterMenuItems>
      <FilterMenuItem value={labelMaterial} depth={0}/>
      <FilterMenuItem value={labelElements} depth={1}/>
      <FilterMenuItem value={labelSymmetry} depth={1}/>
      <FilterMenuItem value={labelMethod} depth={0}/>
      <FilterMenuItem value={labelSimulation} depth={1}/>
      <FilterMenuItem value={labelDFT} depth={2}/>
      <FilterMenuItem value={labelGW} depth={2}/>
      <FilterMenuItem value={labelExperiment} depth={1} disableButton/>
      <FilterMenuItem value={labelEELS} depth={2}/>
      <FilterMenuItem value={labelProperties} depth={0} disableButton/>
      <FilterMenuItem value={labelElectronic} depth={1}/>
      <FilterMenuItem value={labelVibrational} depth={1}/>
      <FilterMenuItem value={labelMechanical} depth={1}/>
      <FilterMenuItem value={labelSpectroscopy} depth={1}/>
      <FilterMenuItem value={labelGeometryOptimization} depth={1}/>
      <FilterMenuItem value={labelELN} depth={0}/>
      <FilterMenuItem value={labelAuthor} depth={0}/>
      <FilterMenuItem value={labelDataset} depth={0}/>
      <FilterMenuItem value={labelAccess} depth={0}/>
      <FilterMenuItem value={labelIDs} depth={0}/>
      <FilterMenuItem value={labelArchive} depth={0}/>
      {resource === 'materials' &&
        <InputCheckbox
          quantity="combine"
          label="Combine results from several entries"
          description="If selected, your filters may be matched from several
          entries that contain the same material. When unchecked, the material
          has to have a single entry that matches all your filters."
          className={styles.combine}
        ></InputCheckbox>
      }
    </FilterMenuItems>
    <FilterSubMenus>
      {loaded && <>
        <FilterSubMenuMaterial value={labelMaterial}/>
        <FilterSubMenuElements value={labelElements} size="large"/>
        <FilterSubMenuSymmetry value={labelSymmetry}/>
        <FilterSubMenuMethod value={labelMethod}/>
        <FilterSubMenuSimulation value={labelSimulation}/>
        <FilterSubMenuDFT value={labelDFT}/>
        <FilterSubMenuGW value={labelGW}/>
        <FilterSubMenuEELS value={labelEELS}/>
        <FilterSubMenuElectronic value={labelElectronic}/>
        <FilterSubMenuVibrational value={labelVibrational}/>
        <FilterSubMenuMechanical value={labelMechanical}/>
        <FilterSubMenuSpectroscopy value={labelSpectroscopy}/>
        <FilterSubMenuGeometryOptimization value={labelGeometryOptimization}/>
        <FilterSubMenuELN value={labelELN}/>
        <FilterSubMenuAuthor value={labelAuthor} size="medium"/>
        <FilterSubMenuDataset value={labelDataset}/>
        <FilterSubMenuAccess value={labelAccess}/>
        <FilterSubMenuIDs value={labelIDs}/>
        <FilterSubMenuArchive value={labelArchive} size="medium"/>
      </>}
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
