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
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import {
  Menu,
  MenuItem,
  FormControl,
  Tooltip,
  InputLabel,
  Select,
  Button
} from '@material-ui/core'
import { atom, atomFamily, useRecoilState, useRecoilCallback, useRecoilValue } from 'recoil'

const presets = {
  new: {
    icon: 'filled',
    iconSize: 'small',
    menu: 'visible',
    inputVariant: 'filled',
    inputSize: 'small',
    aggIndicator: 'on',
    aggCollapse: 'off'
  },
  old: {
    icon: 'plain',
    iconSize: 'medium',
    menu: 'hidden',
    inputVariant: 'outlined',
    inputSize: 'medium',
    aggIndicator: 'off',
    aggCollapse: 'off'
  }
}

/**
 * Item in the GUI menu.
 */
const usePrototypeMenuItemStyles = makeStyles(theme => ({
  root: {},
  form: {
    width: '100%'
  }
}))
const GUIMenuItem = React.memo(({title, tooltip, value, onChange, options}) => {
  const styles = usePrototypeMenuItemStyles()
  return <MenuItem className={styles.root}>
    <FormControl className={styles.form}>
      <Tooltip title={tooltip || ''}>
        <InputLabel>{title}</InputLabel>
      </Tooltip>
      <Select value={value} onChange={(event) => onChange(event.target.value)}>
        {options.map((option) => <MenuItem key={value} value={option}>{option}</MenuItem>)}
      </Select>
    </FormControl>
  </MenuItem>
})
GUIMenuItem.propTypes = {
  title: PropTypes.string,
  tooltip: PropTypes.string,
  value: PropTypes.string,
  onChange: PropTypes.func,
  options: PropTypes.array
}

/**
 * GUI menu that can be used to dynamically control the GUI layout the menu is
 * hidden by default but can be shown by pressing Alt+G.
 */
const useStyles = makeStyles(theme => ({
  root: {
    position: 'fixed',
    bottom: theme.spacing(2),
    right: theme.spacing(2),
    zIndex: 2
  },
  menu: {
    width: theme.spacing(25)
  }
}))
const GUIMenu = React.memo(() => {
  const styles = useStyles()
  const [visible, setVisible] = useState(false)
  const [icon, setIcon] = useRecoilState(guiState('icon'))
  const [iconSize, setIconSize] = useRecoilState(guiState('iconSize'))
  const [menu, setMenu] = useRecoilState(guiState('menu'))
  const [inputVariant, setInputVariant] = useRecoilState(guiState('inputVariant'))
  const [inputSize, setInputSize] = useRecoilState(guiState('inputSize'))
  const [aggIndicator, setAggIndicator] = useRecoilState(guiState('aggIndicator'))
  const [aggCollapse, setAggCollapse] = useRecoilState(guiState('aggCollapse'))
  const preset = useRecoilValue(presetState)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const open = Boolean(anchorEl)
  const handleClick = (event) => {
    setAnchorEl(event.currentTarget)
  }
  const handleClose = () => {
    setAnchorEl(null)
  }

  // When the menu is mounted, add a global key listener that can be used to
  // toggle the menu.
  useEffect(() => {
    function handleHotkey(event) {
      if (event.altKey && event.key === 'g') {
        setAnchorEl(old => old === null ? old : null)
        setVisible(old => !old)
      }
    }
    window.addEventListener('keydown', handleHotkey)
    return () => window.removeEventListener('keydown', handleHotkey)
  }, [])

  const updatePreset = useRecoilCallback(({set}) => (preset) => {
    const presetMap = presets[preset]
    for (const [key, value] of Object.entries(presetMap)) {
      set(guiState(key), value)
    }
    set(presetState, preset)
  }, [])

  return <div className={styles.root} style={{display: visible ? 'block' : 'none'}}>
    <Button variant="contained" color="primary" onClick={handleClick}>
      GUI options
    </Button>
    <Menu
      anchorEl={anchorEl}
      open={open}
      onClose={handleClose}
    >
      <div className={styles.menu}>
        <GUIMenuItem
          title="Preset"
          value={preset}
          onChange={updatePreset}
          tooltip="Switch to a specific preset style."
          options={['new', 'old']}
        />
        <GUIMenuItem
          title="Text input variant"
          value={inputVariant}
          onChange={setInputVariant}
          tooltip="The MUI variant for text input fields."
          options={['outlined', 'filled', 'standard']}
        />
        <GUIMenuItem
          title="Text input size"
          value={inputSize}
          onChange={setInputSize}
          tooltip="The MUI size for text input fields."
          options={['small', 'medium']}
        />
        <GUIMenuItem
          title="Input options loading indicator"
          value={aggIndicator}
          onChange={setAggIndicator}
          tooltip="Whether to show a loading indicator when input field options are populated for the first time."
          options={['on', 'off']}
        />
        <GUIMenuItem
          title="Input options collapse"
          value={aggCollapse}
          onChange={setAggCollapse}
          tooltip="Controls if there is a fixed space reserved for input field options, or whether the space is determined by the amount of items that are shown."
          options={['on', 'off']}
        />
        <GUIMenuItem
          title="Filter settings display"
          value={menu}
          onChange={setMenu}
          tooltip="Controls the display style for filter settings that include e.g. statistics scaling."
          options={['hidden', 'visible']}
        />
        <GUIMenuItem
          title="Statistics icon style"
          value={icon}
          onChange={setIcon}
          tooltip="Controls the style of the icon that is used to add/remove a filter from the statistics grid."
          options={['filled', 'outlined', 'plain']}
        />
        <GUIMenuItem
          title="Statistics icon size"
          value={iconSize}
          onChange={setIconSize}
          tooltip="Controls the size of the icon that is used to add/remove a filter from the statistics grid."
          options={['small', 'medium']}
        />
      </div>
    </Menu>
  </div>
})
const presetDefault = 'new'
export const presetState = atom({
  key: 'preset',
  default: presetDefault
})
export const guiState = atomFamily({
  key: 'guiState',
  default: (name) => presets[presetDefault][name]
})

export default GUIMenu
