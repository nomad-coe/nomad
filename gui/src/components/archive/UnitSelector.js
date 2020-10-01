import React, { useCallback, useState } from 'react'
import { useRecoilState } from 'recoil'
import { makeStyles } from '@material-ui/core/styles'
import {
  Box,
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
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { conversionMap, unitMap, unitSystems } from '../../units'

/**
 * Component that wraps it's children in a container that can be 'floated',
 * i.e. displayed on an html element that is positioned relative to the
 * viewport and is above all other elements.
 */
export function UnitSelector({className, classes, unitsState, onUnitChange, onSystemChange}) {
  // States
  const [canSelect, setCanSelect] = useState(true)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const open = Boolean(anchorEl)
  const [units, setUnits] = useRecoilState(unitsState)

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      menuItem: {
        width: '10rem'
      },
      systems: {
        margin: theme.spacing(2),
        marginTop: theme.spacing(1)
      }
    }
  })
  const style = useStyles(classes)

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
    console.log(changes)
    setUnits({...units, ...changes})
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Ordered list of controllable units. It may be smaller than the full list of
  // units.
  const unitNames = ['energy', 'length', 'force', 'mass', 'time', 'temperature']
  const systemNames = ['SI', 'AU']

  return (
    <Box className={clsx(style.root, className)}>
      <Button
        aria-controls="customized-menu"
        aria-haspopup="true"
        variant="outlined"
        color="primary"
        onClick={openMenu}
      >
        Select units
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
          classes={{root: style.systems}}
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
                classes={{root: style.menuItem}}
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
    </Box>
  )
}

UnitSelector.propTypes = {
  /**
   * CSS class for the root element.
   */
  className: PropTypes.string,
  /**
   * CSS classes for this component.
   */
  classes: PropTypes.object,
  /**
   * Recoil atom containing the unit configuration that this component will
   * attach to.
   */
  unitsState: PropTypes.object,
  /**
   * Callback for unit selection.
   */
  onUnitChange: PropTypes.func,
  /**
   * Callback for unit system selection.
   */
  onSystemChange: PropTypes.func
}
UnitSelector.defaultProps = {
  float: false
}
