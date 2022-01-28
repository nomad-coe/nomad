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
import React, { useCallback, useState } from 'react'
import { useRecoilState } from 'recoil'
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
import { unitsState } from '../units'
import { conversionMap, unitMap, unitSystems } from '../unitsData'

/**
 * Component that wraps it's children in a container that can be 'floated',
 * i.e. displayed on an html element that is positioned relative to the
 * viewport and is above all other elements.
 */
const useStyles = makeStyles((theme) => {
  return {
    root: {
    },
    menuItem: {
      width: '10rem'
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
  const [canSelect, setCanSelect] = useState(true)
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
  const handleSystemChange = useCallback((event) => {
    const systemName = event.target.value
    let changes = {system: systemName}
    if (systemName === 'custom') {
      setCanSelect(true)
    } else {
      setCanSelect(false)
      const system = unitSystems[systemName]
      changes = {...changes, ...system.units}
    }
    setUnits({...units, ...changes})
    if (onSystemChange) {
      onSystemChange(event)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])
  const handleUnitChange = useCallback(event => {
    const changes = {[event.target.name]: event.target.value}
    if (onUnitChange) {
      onUnitChange(event)
    }
    setUnits({...units, ...changes})
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Ordered list of controllable units. It may be smaller than the full list of
  // units.
  const unitNames = ['energy', 'length', 'force', 'mass', 'time', 'temperature', 'pressure', 'angle']
  const systemNames = ['SI', 'AU']

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
        <RadioGroup aria-label="gender" name="gender1" value={units.system} onChange={handleSystemChange}>
          <Tooltip title="Custom units">
            <FormControlLabel value="custom" control={<Radio />} label="Custom" />
          </Tooltip>
          {systemNames.map((systemName) => {
            const system = unitSystems[systemName]
            return <Tooltip key={systemName} title={system.description}>
              <FormControlLabel value={systemName} control={<Radio />} label={system.label} />
            </Tooltip>
          })}
        </RadioGroup>
      </FormControl>
      {unitNames.map((dimension) => {
        const unitList = conversionMap[dimension].units
        return <MenuItem
          key={dimension}
        >
          <FormControl disabled={!canSelect}>
            <InputLabel id="demo-simple-select-label">{dimension}</InputLabel>
            <Select
              classes={{root: styles.menuItem}}
              labelId="demo-simple-select-label"
              id="demo-simple-select"
              name={dimension}
              value={units[dimension]}
              onChange={handleUnitChange}
            >
              {unitList.map((unit) => {
                const unitLabel = unitMap[unit].label
                return <MenuItem key={unit} value={unit}>{unitLabel}</MenuItem>
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
UnitSelector.defaultProps = {
  float: false
}

export default UnitSelector
