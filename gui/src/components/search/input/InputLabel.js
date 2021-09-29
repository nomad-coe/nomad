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
import React, { useMemo, useCallback, useState } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import {
  Tooltip,
  Typography,
  Radio,
  FormControl,
  FormLabel,
  FormControlLabel,
  RadioGroup,
  Menu
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import AddIcon from '@material-ui/icons/Add'
import RemoveIcon from '@material-ui/icons/Remove'
import DragHandleIcon from '@material-ui/icons/DragHandle'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import { useStatisticState } from '../SearchContext'
import { Actions, Action } from '../../Actions'

/**
 * The quantity label and actions shown by all filter components.
 */
const useStaticStyles = makeStyles(theme => ({
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
  label: {
    textTransform: 'capitalize',
    fontSize: '0.9rem',
    color: '#383838'
  }
}))
const scales = {
  'linear': 1,
  '1/2': 0.5,
  '1/4': 0.25,
  '1/8': 0.125
}

const FilterLabel = React.memo(({
  quantity,
  label,
  underscores,
  description,
  disableStatistics,
  scale,
  onChangeScale,
  draggable,
  className,
  classes
}) => {
  const styles = useStaticStyles({classes: classes})
  const [statistic, setStatistic] = useStatisticState(quantity)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const isSettingsOpen = Boolean(anchorEl)
  const [isTooltipOpen, setIsTooltipOpen] = useState(false)
  const [isDragging, setIsDragging] = useState(false)

  // Callbacks
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])
  const openTooltip = useCallback((event) => {
    !isDragging && setIsTooltipOpen(true)
  }, [isDragging])
  const closeTooltip = useCallback(() => {
    setIsTooltipOpen(false)
  }, [])
  const handleMouseDown = useCallback((event) => {
    setIsDragging(true)
    setIsTooltipOpen(false)
  }, [])
  const handleMouseUp = useCallback(() => {
    setIsDragging(false)
  }, [])

  const tooltipProps = useMemo(() => ({
    open: isTooltipOpen,
    onClose: closeTooltip,
    onOpen: openTooltip
  }), [closeTooltip, openTooltip, isTooltipOpen])

  // Remove underscores from name
  const finalLabel = useMemo(() => {
    let finalLabel = label || quantity
    return !underscores ? finalLabel.replace(/_/g, ' ') : finalLabel
  }, [label, quantity, underscores])

  // The tooltip needs to be controlled: otherwise it won't close as we open the
  // select menu
  return <Actions
    className={clsx(className, styles.root)}
    header={
      <Tooltip title={description || ''} placement="bottom">
        <Typography
          className={styles.label}
          variant="button"
        >
          {finalLabel}
        </Typography>
      </Tooltip>
    }
  >
    {!disableStatistics && <>
      <Action
        tooltip={statistic ? 'Remove the filter statistics from the results panel.' : 'Show filter statistics in the results panel.'}
        disabled={disableStatistics}
        onClick={() => setStatistic(old => !old)}
      >
        {statistic ? <RemoveIcon/> : <AddIcon/>}
      </Action>
      {draggable && <Action
        onMouseDown={handleMouseDown}
        onMouseUp={handleMouseUp}
        className="dragHandle"
        tooltip="Drag to change location"
        TooltipProps={tooltipProps}
      >
        <DragHandleIcon/>
      </Action>
      }
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
            onChange={(event, value) => onChangeScale(Number(value))}
          >
            {Object.entries(scales).map(([key, value]) =>
              <FormControlLabel key={key} value={value} label={key} control={<Radio/>} />
            )}
          </RadioGroup>
        </FormControl>
      </Menu>
    </>}
  </Actions>
})

FilterLabel.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  underscores: PropTypes.bool,
  disableStatistics: PropTypes.bool,
  disableAggSize: PropTypes.bool,
  scale: PropTypes.oneOf(Object.values(scales)),
  onChangeScale: PropTypes.func,
  draggable: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object
}

FilterLabel.defaultProps = {
  underscores: false,
  disableStatistics: false,
  scale: 1
}

export default FilterLabel
