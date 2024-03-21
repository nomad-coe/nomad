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
import React, { useCallback, useState, useMemo } from 'react'
import {
  Box, Button, Menu, FormLabel, makeStyles, Typography
} from '@material-ui/core'
import SettingsIcon from '@material-ui/icons/Settings'
import ReplayIcon from '@material-ui/icons/Replay'
import PropTypes from 'prop-types'
import { HelpButton } from '../Help'
import { InputText } from '../search/input/InputText'
import UnitDimensionSelect from './UnitDimensionSelect'
import UnitSystemSelect from './UnitSystemSelect'
import { useUnitContext } from './UnitContext'
import { Action, ActionHeader, Actions } from '../Actions'

/**
 * Menu for controlling all units in the current unit context.
 */
const useStyles = makeStyles(theme => ({
  // MUI will automatically add a padding whe scroll bar is visible. This is
  // disabled here because the contents already have a sufficient margin.
  list: {
    paddingRight: '0 !important',
    width: '100% !important'
  },
  option: {
    paddingTop: '6px',
    paddingBottom: '6px'
  }
}))
const UnitMenu = React.memo(({
  className,
  onUnitChange,
  onSystemChange
}) => {
  const {units, dimensionMap, reset} = useUnitContext()
  const [anchorEl, setAnchorEl] = useState(null)
  const open = Boolean(anchorEl)
  const styles = useStyles()

  const dimensionOptions = useMemo(() => {
    return Object.keys(units)
      .filter((name) => name !== 'dimensionless')
      .sort()
  }, [units])
  const [dimension, setDimension] = useState(dimensionOptions[0])
  const [dimensionInput, setDimensionInput] = useState(dimensionMap[dimensionOptions[0]].label)

  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

  const handleChange = useCallback((value) => {
    setDimensionInput(value)
  }, [])

  const handleBlur = useCallback(() => {
    const cleanValue = dimensionInput.trim().toLowerCase()
    const dim = dimensionMap[cleanValue]
    if (dim) {
      setDimensionInput(dim.label)
      setDimension(cleanValue)
    } else {
      setDimensionInput(dimensionMap[dimension].label)
    }
  }, [dimension, dimensionInput, dimensionMap])

  const handleSelect = useCallback((value) => {
    setDimension(value)
    setDimensionInput(dimensionMap[value].label)
  }, [dimensionMap])

  return <>
    <Button
      aria-controls="customized-menu"
      aria-haspopup="true"
      variant="text"
      color="primary"
      onClick={openMenu}
      className={className}
      startIcon={<SettingsIcon/>}
    >
      Units
    </Button>
    <Menu
      variant="menu"
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      keepMounted
      open={open}
      onClose={closeMenu}
      classes={{list: styles.list}}
    >
      <Box px={2.5} py={1} width="20rem">
        <Actions>
          <ActionHeader>
            <Typography variant='h6' fontSize='0.9rem'>
              Unit Settings
            </Typography>
          </ActionHeader>
          <Action
            tooltip="About units"
            ButtonComponent={HelpButton}
            ButtonProps={{
              title: "About units",
              IconProps: {fontSize: 'small'},
              content: `
                With these settings you can change in which units numerical data
                is **displayed** in the browser. This display unit can be
                different from the unit that is used when the data is stored.
                Note that it is possible to define a fixed display unit in the
                metainfo which will overrule the settings made through this
                menu.

                Each NOMAD installation comes with a default set of unit
                systems, which can be modified in the \`nomad.yaml\`
                configuration file. You can here choose which of these unit
                systems to use.

                Each unit system contains information about the exact unit that
                should be used for each dimension. Many of the commonly used
                units are available for selection, and you may use the SI
                prefixes on any of them. Note that in some unit systems certain
                dimensions are locked, which means that you cannot change them.
                `
            }}
          >
          </Action>
          <Action tooltip="Reset unit settings" onClick={reset}>
            <ReplayIcon fontSize="small"/>
          </Action>
        </Actions>
        <Box mt={1} />
        {/* <FormControl variant="filled" fullWidth disabled={!scope}>
          <InputLabel>Scope</InputLabel>
          <Select
            value={scope || 'Global'}
            onChange={(event) => setScope(event.target.value)}
          >
            <MenuItem value={'Global'}>Global</MenuItem>
            <MenuItem value={'Schema'}>Schema</MenuItem>
          </Select>
        </FormControl> */}
        <Box mt={1} />
        <FormLabel component="legend">Select unit system</FormLabel>
        <UnitSystemSelect onChange={onSystemChange}/>
        <Box mt={1} />
        <FormLabel component="legend">Select dimension and unit</FormLabel>
        <Box mt={1} />
        <InputText
          value={dimensionInput}
          suggestions={dimensionOptions}
          disableClearable
          suggestAllOnFocus
          showOpenSuggestions
          onChange={handleChange}
          onSelect={handleSelect}
          onBlur={handleBlur}
          renderOption={(option) => <Typography className={styles.option}>
            {dimensionMap[option].label}
          </Typography>}
          getOptionLabel={(option) => option}
          TextFieldProps={{label: 'Dimension'}}
        />
        <Box mt={1} />
        <UnitDimensionSelect
          onChange={onUnitChange}
          dimension={dimension}
          label="Unit"
        />
      </Box>
    </Menu>
  </>
})

UnitMenu.propTypes = {
  onUnitChange: PropTypes.func,
  onSystemChange: PropTypes.func,
  className: PropTypes.string
}

export default UnitMenu
