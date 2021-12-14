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
import { makeStyles } from '@material-ui/core/styles'
import {
  Checkbox,
  Radio,
  Tooltip,
  Typography,
  FormControlLabel
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useSearchContext } from '../SearchContext'
import StatisticsBar from '../statistics/StatisticsBar'

/**
 * Represents a selectable item for a filter value.
*/
export const inputItemHeight = 34
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: inputItemHeight,
    position: 'relative'
  },
  bar: {
    position: 'absolute',
    top: theme.spacing(0.8),
    left: theme.spacing(3.2),
    right: 0,
    bottom: theme.spacing(0.8)
  },
  controlLabel: {
    position: 'absolute',
    top: 0,
    left: 0,
    bottom: 0,
    right: '4rem'
  },
  label: {
    width: '100%'
  }
}))
const InputItem = React.memo(({
  value,
  label,
  selected,
  onChange,
  disabled,
  tooltip,
  variant,
  total,
  count,
  scale,
  disableStatistics,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = useStyles(classes)
  const { useIsStatisticsEnabled } = useSearchContext()
  const isStatisticsEnabled = useIsStatisticsEnabled()

  const handleChange = useCallback((event, itemValue) => {
    if (!disabled && onChange) onChange(value, itemValue)
  }, [value, disabled, onChange])

  let Control
  if (variant === 'radio') {
    Control = Radio
  } else if (variant === 'checkbox') {
    Control = Checkbox
  }

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    {(isStatisticsEnabled && !disableStatistics) && <StatisticsBar
      className={styles.bar}
      max={total}
      value={count}
      scale={scale}
      selected={selected}
      disabled={disabled}
    />}
    <FormControlLabel
      className={styles.controlLabel}
      classes={{label: styles.label}}
      disabled={disabled}
      control={<Control
        checked={selected}
        size="medium"
        color="primary"
        onChange={handleChange}
        name={value}
      />}
      label={
        <Tooltip
          placement="right"
          enterDelay={200}
          title={tooltip || ''}
        >
          <Typography className={styles.label} noWrap>{label || value}</Typography>
        </Tooltip>
      }
    />
  </div>
})

InputItem.propTypes = {
  value: PropTypes.string, // The actual value
  label: PropTypes.string, // The name to show
  selected: PropTypes.bool, // Whether the option is selected or not
  onChange: PropTypes.func, // Callback when selecting
  disabled: PropTypes.bool, // Whether the option should be disabled
  tooltip: PropTypes.string, // Tooltip that is shown for label
  variant: PropTypes.oneOf(['radio', 'checkbox']), // The type of item to display
  total: PropTypes.number, // Total number for statistics
  count: PropTypes.number, // Count of these values for statistics
  scale: PropTypes.number, // Scaling of the statistics
  disableStatistics: PropTypes.bool, // Use to disable statistics for this item
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default InputItem
