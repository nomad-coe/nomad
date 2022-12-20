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
  Menu
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useRecoilValue } from 'recoil'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import FilterTitle from '../FilterTitle'
import { Actions, ActionHeader, Action, ActionSelect } from '../../Actions'
import WidgetToggle from '../widgets/WidgetToggle'
import { scales } from '../../plotting/common'
import { guiState } from '../../GUIMenu'
import { useBoolState } from '../../../hooks'

/**
 * The quantity label and actions shown by all filter components.
 */
const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(0.5),
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
    flexGrow: 1,
    minWidth: 0
  },
  spacer: {
    flexGrow: 1,
    minWidth: 0
  }
}))
const InputHeader = React.memo(({
  quantity,
  label,
  description,
  disableWidget,
  disableStatistics,
  scale,
  onChangeScale,
  actions,
  actionsAlign,
  className,
  classes
}) => {
  const styles = useStyles({classes: classes})
  const [anchorEl, setAnchorEl] = React.useState(null)
  const isSettingsOpen = Boolean(anchorEl)
  const [isDragging, setDragging, setNotDragging] = useBoolState(false)
  const [isTooltipOpen, openTooltip, closeTooltip] = useBoolState(false)
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
            onChange={onChangeScale ? (event, value) => onChangeScale(value) : undefined}
          >
            {Object.entries(scales).map(([key, value]) =>
              <FormControlLabel key={key} value={key} label={key} control={<Radio/>} />
            )}
          </RadioGroup>
        </FormControl>
      </Menu>
    </>
    : <ActionSelect
      value={scale}
      options={Object.keys(scales)}
      tooltip="Statistics scaling"
      onChange={onChangeScale}
    />

  return <Actions className={clsx(styles.root, className)}>
    <ActionHeader disableSpacer>
      <div
        className={clsx(styles.row)}
        onMouseDown={handleMouseDown}
        onMouseUp={handleMouseUp}
      >
        <FilterTitle
          quantity={quantity}
          label={label}
          description={description}
          TooltipProps={tooltipProps}
        />
        <div className={styles.spacer}/>
      </div>
    </ActionHeader>
    {actionsAlign === 'left' && actions}
    {!disableStatistics && menuComp}
    {!disableWidget && <WidgetToggle quantity={quantity} disabled={disableWidget} />}
    {actionsAlign === 'right' && actions}
  </Actions>
})

InputHeader.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  disableWidget: PropTypes.bool,
  disableStatistics: PropTypes.bool,
  scale: PropTypes.string,
  onChangeScale: PropTypes.func,
  variant: PropTypes.string,
  actions: PropTypes.node,
  actionsAlign: PropTypes.oneOf(['left', 'right']),
  className: PropTypes.string,
  classes: PropTypes.object
}

InputHeader.defaultProps = {
  underscores: false,
  disableWidget: false,
  disableStatistics: false,
  actionsAlign: 'left',
  scale: 'linear'
}

export default InputHeader
