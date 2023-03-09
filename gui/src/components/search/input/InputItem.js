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
import StatisticsBar from '../widgets/StatisticsBar'

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
  controlLabel: {
    width: '100%',
    height: '100%',
    margin: 0
  },
  labelContainer: {
    height: '100%',
    width: '100%'
  },
  labelPlacementStart: {
    margin: 0
  },
  container: {
    position: 'relative',
    height: '100%',
    width: '100%'
  },
  label: {
    position: 'absolute',
    top: 0,
    left: theme.spacing(0.5),
    right: 0,
    bottom: 0,
    lineHeight: '100%',
    textAlign: 'center',
    display: 'flex',
    alignItems: 'center'
  },
  bar: {
    position: 'absolute',
    top: theme.spacing(0.8),
    left: 0,
    right: 0,
    bottom: theme.spacing(0.8)
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
  max,
  count,
  scale,
  disableStatistics,
  labelPlacement,
  disableLabelClick,
  disableSelect,
  className,
  classes,
  'data-testid': testID,
  ...moreProps
}) => {
  const styles = useStyles(classes)

  const finalLabel = label || `${value}`
  const handleChange = useCallback((event) => {
    if (!disabled && onChange) onChange(event, value, !selected)
  }, [value, disabled, onChange, selected])

  let Control
  if (variant === 'radio') {
    Control = Radio
  } else if (variant === 'checkbox') {
    Control = Checkbox
  }

  // Component that contains the label and the statistics
  const labelComponent = <div className={styles.container}>
    {!disableStatistics && <StatisticsBar
      className={styles.bar}
      max={max}
      value={count}
      scale={scale}
      selected={selected}
      disabled={disabled}
    />}
    <div className={styles.label} style={!disableStatistics ? {right: '4rem'} : undefined}>
      <Tooltip
        placement="right"
        enterDelay={200}
        title={tooltip || ''}
      >
        <Typography noWrap>{finalLabel}</Typography>
      </Tooltip>
    </div>
  </div>

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    {disableSelect
      ? labelComponent
      : <FormControlLabel
        style={{pointerEvents: disableLabelClick ? 'none' : undefined}}
        className={styles.controlLabel}
        classes={{
          label: styles.labelContainer,
          labelPlacementStart: styles.labelPlacementStart
        }}
        labelPlacement={labelPlacement}
        disabled={disabled}
        control={<Control
          checked={selected}
          size="medium"
          color="secondary"
          onClick={handleChange}
          name={finalLabel}
          style={{
            padding: 6,
            pointerEvents: (disableLabelClick ? 'auto' : undefined),
            marginRight: (labelPlacement === 'start' ? -8 : undefined),
            marginLeft: (labelPlacement === 'end' ? -8 : undefined)
          }}
          {...moreProps}
        />}
        label={labelComponent}
      />
    }
  </div>
})

InputItem.propTypes = {
  value: PropTypes.any, // The actual value
  label: PropTypes.string, // The name to show
  selected: PropTypes.bool, // Whether the option is selected or not
  onChange: PropTypes.func, // Callback when selecting
  disabled: PropTypes.bool, // Whether the option should be disabled
  tooltip: PropTypes.string, // Tooltip that is shown for label
  variant: PropTypes.oneOf(['radio', 'checkbox']), // The type of item to display
  max: PropTypes.number, // Maximum for statistics
  count: PropTypes.number, // Count of these values for statistics
  scale: PropTypes.string, // Scaling of the statistics
  disableStatistics: PropTypes.bool, // Use to disable statistics for this item
  labelPlacement: PropTypes.oneOf(['start', 'end']), // Controls the label placement
  disableLabelClick: PropTypes.bool, // Whether clicking the label is enabled
  disableSelect: PropTypes.bool, // Whether the checkbox is shown
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

InputItem.defaultProps = {
  labelPlacement: 'end'
}

export default InputItem
