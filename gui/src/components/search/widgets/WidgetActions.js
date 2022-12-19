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
import React from 'react'
import {
  MenuItem,
  Select,
  Checkbox,
  Tooltip,
  FormControlLabel
} from '@material-ui/core'
import { isArray } from 'lodash'
import PropTypes from 'prop-types'
import { Action } from '../../Actions'
import { useBoolState } from '../../../hooks'

/**
 * Dropdown menu for controlling widget state
 */
export const WidgetActionSelect = React.memo(({value, options, tooltip, onChange}) => {
  const [isStatsTooltipOpen, openStatsTooltip, closeStatsTooltip] = useBoolState(false)
  const items = isArray(options)
    ? Object.fromEntries(options.map(x => [x, x]))
    : options
  return <Action
    TooltipProps={{
      title: tooltip || "",
      open: isStatsTooltipOpen,
      disableHoverListener: true
  }}>
    <Select
      value={value}
      onMouseEnter={openStatsTooltip}
      onMouseLeave={closeStatsTooltip}
      onOpen={closeStatsTooltip}
      onChange={(event) => onChange && onChange(event.target.value)}
    >
      {Object.entries(items).map(([key, value]) =>
        <MenuItem key={key} value={value}>{key}</MenuItem>
      )}
    </Select>
  </Action>
})

WidgetActionSelect.propTypes = {
  value: PropTypes.string,
  options: PropTypes.oneOfType([PropTypes.object, PropTypes.arrayOf(PropTypes.string)]),
  tooltip: PropTypes.string,
  onChange: PropTypes.func
}

/**
 * Checkbox for controlling widget state
 */
export const WidgetActionCheckbox = React.memo(({value, label, tooltip, onChange}) => {
  return <Tooltip title={tooltip}>
      <FormControlLabel
        control={<Checkbox
          checked={value}
          onChange={(event, value) => onChange && onChange(value)}
          size="small"
        />}
        label={label}
      />
    </Tooltip>
})

WidgetActionCheckbox.propTypes = {
  value: PropTypes.bool,
  label: PropTypes.string,
  tooltip: PropTypes.string,
  onChange: PropTypes.func
}
