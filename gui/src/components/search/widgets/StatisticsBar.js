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
import React, { useMemo } from 'react'
import clsx from 'clsx'
import { Typography, Tooltip } from '@material-ui/core/'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import { approxInteger } from '../../../utils'
import { getScaler } from '../../plotting/common'

/**
 * A rectangular bar displaying the relative occurence of a specific value. Uses
 * an animated CSS transform to react to changes in the bar size. Notice that
 * using a transform is much more performant than animating the width.
 */
const useStyles = makeStyles(theme => ({
  root: {
  },
  container: {
    position: 'relative',
    width: '100%',
    height: '100%'
  },
  rectangle: {
    backgroundColor: theme.palette.primary.light,
    width: '100%',
    height: '100%',
    '-webkit-transform': 'none',
    transform: 'none',
    transition: 'transform 250ms',
    transformOrigin: 'bottom left',
    willChange: 'transform'
  },
  value: {
    position: 'absolute',
    right: theme.spacing(1),
    top: 0,
    bottom: 0
  }
}))
const StatisticsBar = React.memo(({
  max,
  value,
  scale,
  selected,
  vertical,
  disabled,
  disableValue,
  tooltip,
  onClick,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = useStyles(classes)
  const theme = useTheme()

  // Calculate the approximated count and the final scaled value
  const scaler = useMemo(() => scale ? getScaler(scale) : (value) => value, [scale])
  const finalCount = useMemo(() => approxInteger(value || 0), [value])
  const finalScale = useMemo(() => scaler(value / max) || 0, [value, max, scaler])

  return <div onClick={onClick} className={clsx(className, styles.root)} data-testid={testID}>
    <Tooltip placement="bottom" enterDelay={0} title={tooltip || ''}>
      <div className={styles.container}>
        <div className={styles.rectangle} style={{
          transform: vertical ? `scaleY(${finalScale})` : `scaleX(${finalScale})`,
          backgroundColor: disabled ? theme.palette.action.disabledBackground : (selected ? theme.palette.secondary.veryLight : theme.palette.primary.veryLight)
        }}></div>
        {!disableValue && <Typography className={styles.value} style={{color: disabled ? theme.palette.text.disabled : undefined}}>{finalCount}</Typography>}
      </div>
    </Tooltip>
  </div>
})

StatisticsBar.propTypes = {
  max: PropTypes.number,
  value: PropTypes.number,
  scale: PropTypes.string,
  selected: PropTypes.bool,
  disabled: PropTypes.bool,
  disableValue: PropTypes.bool,
  tooltip: PropTypes.string,
  vertical: PropTypes.bool,
  onClick: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default StatisticsBar
