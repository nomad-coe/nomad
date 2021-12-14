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
import React, { useMemo, useCallback } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import {
  Radio,
  FormControl,
  FormLabel,
  FormControlLabel,
  RadioGroup,
  Menu,
  MenuItem,
  Select
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useRecoilValue } from 'recoil'
import MUIAddIcon from '@material-ui/icons/Add'
import AddCircleOutlineIcon from '@material-ui/icons/AddCircleOutline'
import CloseIcon from '@material-ui/icons/Close'
import HighlightOffIcon from '@material-ui/icons/HighlightOff'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import AddCircleIcon from '@material-ui/icons/AddCircle'
import CancelIcon from '@material-ui/icons/Cancel'
import InputTitle from './InputTitle'
import { useSearchContext } from '../SearchContext'
import { Actions, ActionHeader, Action } from '../../Actions'
import { guiState } from '../../GUIMenu'
import { useBoolState } from '../../../hooks'

/**
 * The quantity label and actions shown by all filter components.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(0.5),
    marginTop: theme.spacing(-0.5),
    height: '2.5rem',
    width: '100%'
  },
  menuItem: {
    width: '10rem',
    margin: theme.spacing(2),
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1)
  },
  row: {
    display: 'flex',
    height: '100%',
    alignItems: 'center',
    flexGrow: 1
  },
  spacer: {
    flexGrow: 1
  }
}))
const scales = {
  'linear': 1,
  '1/2': 0.5,
  '1/4': 0.25,
  '1/8': 0.125
}

const InputHeader = React.memo(({
  quantity,
  label,
  description,
  disableStatistics,
  disableScale,
  scale,
  onChangeScale,
  draggable,
  className,
  classes
}) => {
  const { useStatisticState } = useSearchContext()
  const styles = useStaticStyles({classes: classes})
  const [statistic, setStatistic] = useStatisticState(quantity)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const isSettingsOpen = Boolean(anchorEl)
  const [isStatsTooltipOpen, openStatsTooltip, closeStatsTooltip] = useBoolState(false)
  const [isDragging, setDragging, setNotDragging] = useBoolState(false)
  const [isTooltipOpen, openTooltip, closeTooltip] = useBoolState(false)
  const align = useRecoilValue(guiState('align'))
  const icon = useRecoilValue(guiState('icon'))
  const iconSize = useRecoilValue(guiState('iconSize'))
  const menu = useRecoilValue(guiState('menu'))

  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])
  const handleMouseDown = useCallback((event) => {
    setDragging()
    closeTooltip()
  }, [closeTooltip, setDragging])
  const handleMouseUp = useCallback(() => {
    setNotDragging()
  }, [setNotDragging])
  const tooltipProps = useMemo(() => ({
    open: isTooltipOpen,
    onClose: closeTooltip,
    onOpen: () => !isDragging && openTooltip()
  }), [isTooltipOpen, closeTooltip, isDragging, openTooltip])

  let RemoveIcon, AddIcon
  if (icon === 'plain') {
    RemoveIcon = CloseIcon
    AddIcon = MUIAddIcon
  } else if (icon === 'outlined') {
    RemoveIcon = HighlightOffIcon
    AddIcon = AddCircleOutlineIcon
  } else if (icon === 'filled') {
    RemoveIcon = CancelIcon
    AddIcon = AddCircleIcon
  }

  const control = <Action
    tooltip={statistic ? 'Remove the filter statistics from the results panel.' : 'Show filter statistics in the results panel.'}
    disabled={disableStatistics}
    onClick={() => setStatistic(old => !old)}
  >
    {statistic ? <RemoveIcon fontSize={iconSize}/> : <AddIcon fontSize={iconSize}/>}
  </Action>

  // Dedice the menu component. The tooltip of scale dropdown needs to be
  // controlled, otherwise it won't close as we open the select menu.
  const menuComp = menu === 'hidden'
    ? <>
      <Action
        tooltip="Options"
        onClick={openMenu}
      >
        <MoreVertIcon/>
      </Action>
      <Menu
        anchorEl={anchorEl}
        open={isSettingsOpen}
        onClose={closeMenu}
        getContentAnchorEl={null}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
        transformOrigin={{ vertical: 'top', horizontal: 'right' }}
        keepMounted
      >
        <FormControl
          className={styles.menuItem}
          component="fieldset"
        >
          <FormLabel component="legend">Statistics scaling</FormLabel>
          <RadioGroup
            value={scale}
            onChange={onChangeScale ? (event, value) => onChangeScale(Number(value)) : undefined}
          >
            {Object.entries(scales).map(([key, value]) =>
              <FormControlLabel key={key} value={value} label={key} control={<Radio/>} />
            )}
          </RadioGroup>
        </FormControl>
      </Menu>
    </>
    : <Action TooltipProps={{title: 'Statistics scaling', open: isStatsTooltipOpen, disableHoverListener: true}}>
      <Select
        value={scale}
        onMouseEnter={openStatsTooltip}
        onMouseLeave={closeStatsTooltip}
        onOpen={closeStatsTooltip}
        onChange={onChangeScale ? (event) => onChangeScale(Number(event.target.value)) : undefined}
      >
        {Object.entries(scales).map(([key, value]) =>
          <MenuItem key={key} value={value}>{key}</MenuItem>
        )}
      </Select>
    </Action>

  return <Actions className={clsx(className, styles.root)}>
    {align === 'left' && !disableStatistics && control}
    <ActionHeader disableSpacer>
      <div
        className={clsx(styles.row, draggable ? 'dragHandle' : undefined)}
        style={{cursor: draggable ? (isDragging ? 'grabbing' : 'grab') : undefined}}
        onMouseDown={handleMouseDown}
        onMouseUp={handleMouseUp}
      >
        <InputTitle
          quantity={quantity}
          label={label}
          description={description}
          TooltipProps={tooltipProps}
        />
        <div className={styles.spacer}/>
      </div>
    </ActionHeader>
    {!disableScale && menuComp}
    {align === 'right' && !disableStatistics && control}
  </Actions>
})

InputHeader.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  disableStatistics: PropTypes.bool,
  disableScale: PropTypes.bool,
  scale: PropTypes.oneOf(Object.values(scales)),
  onChangeScale: PropTypes.func,
  draggable: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object
}

InputHeader.defaultProps = {
  underscores: false,
  disableStatistics: false,
  disableScale: false,
  scale: 1
}

export default InputHeader
