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
import React, { useCallback } from 'react'
import { useRecoilState } from 'recoil'
import { isNil } from 'lodash'
import { makeStyles } from '@material-ui/core/styles'
import {
  Button,
  Menu,
  MenuItem,
  FormControl,
  InputLabel,
  Select,
  FormLabel,
  FormControlLabel,
  RadioGroup,
  Radio,
  Tooltip
} from '@material-ui/core'
import SettingsIcon from '@material-ui/icons/Settings'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { unitMap, unitSystems, unitsState, dimensionMap } from '../units'

/**
 * Unit selection menu with dropdowns for each dimension and presets for
 * different unit systems.
 */
const useStyles = makeStyles((theme) => {
  return {
    root: {
    },
    menuItem: {
      width: '15rem'
    },
    systems: {
      margin: theme.spacing(2),
      marginTop: theme.spacing(1),
      marginBottom: theme.spacing(1)
    }
  }
})
const UnitSelector = React.memo(({
  className,
  classes,
  onUnitChange,
  onSystemChange
}) => {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const open = Boolean(anchorEl)
  const [units, setUnits] = useRecoilState(unitsState)
  const styles = useStyles({classes: classes})

  // Callbacks
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

  // Used to handle unit system change.
  const handleSystemChange = useCallback((event) => {
    setUnits(unitSystems[event.target.value])
    onSystemChange && onSystemChange(event)
  }, [onSystemChange, setUnits])

  // Used to handle unit change for a specific dimensionality. The changes are
  // stored for each system separately.
  const handleUnitChange = useCallback(event => {
    const dimension = event.target.name
    const unit = event.target.value
    setUnits(old => {
      const newSystem = {
        ...old,
        units: {
          ...old.units,
          ...{[dimension]: {...old.units[dimension], name: unit}}
        }
      }
      unitSystems[old.label] = newSystem
      return unitSystems[old.label]
    })
    onUnitChange && onUnitChange(event)
  }, [onUnitChange, setUnits])

  // Ordered list of controllable units. The 'dimensionless' unit cannot be
  // changed.
  const dimensions = Object.entries(dimensionMap)
    .filter(([dimension, info]) => dimension !== 'dimensionless')

  return <>
    <Button
      aria-controls="customized-menu"
      aria-haspopup="true"
      variant="text"
      color="primary"
      onClick={openMenu}
      className={clsx(styles.root, className)}
      startIcon={<SettingsIcon/>}
    >
      Units
    </Button>
    <Menu
      id="select-unit"
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      keepMounted
      open={open}
      onClose={closeMenu}
    >
      <FormControl
        component="fieldset"
        classes={{root: styles.systems}}
      >
        <FormLabel component="legend">Unit system</FormLabel>
        <RadioGroup aria-label="gender" name="gender1" value={units.label} onChange={handleSystemChange}>
          {Object.values(unitSystems).map(system => {
            return <Tooltip key={system.label} title={system.description}>
              <FormControlLabel value={system.label} control={<Radio />} label={system.label} />
            </Tooltip>
          })}
        </RadioGroup>
      </FormControl>
      {dimensions.map(([dimension, unitInfo]) => {
        const unitDef = units.units[dimension]
        if (isNil(unitDef)) {
          return
        }
        const selectedUnit = unitDef.name
        const disabled = unitDef.fixed
        return <MenuItem
          key={dimension}
        >
          <FormControl disabled={disabled}>
            <InputLabel id="demo-simple-select-label">{unitInfo.label}</InputLabel>
            <Select
              classes={{root: styles.menuItem}}
              labelId="demo-simple-select-label"
              id="demo-simple-select"
              name={dimension}
              value={selectedUnit}
              onChange={handleUnitChange}
            >
              {unitInfo.units.map((unit) => {
                const unitLabel = unitMap[unit].label
                const unitAbbreviation = unitMap[unit].abbreviation
                return <MenuItem key={unit} value={unit}>{`${unitLabel} (${unitAbbreviation})`}</MenuItem>
              })}
            </Select>
          </FormControl>
        </MenuItem>
      })}
    </Menu>
  </>
})

UnitSelector.propTypes = {
  /**
   * Callback for unit selection.
   */
  onUnitChange: PropTypes.func,
  /**
   * Callback for unit system selection.
   */
  onSystemChange: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object
}

export default UnitSelector
