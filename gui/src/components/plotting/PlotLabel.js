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
import clsx from 'clsx'
import { Typography } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import PlotTick from './PlotTick'

/**
 * Label used in plots.
 */
const usePlotLabelStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center'
  },
  horizontal: {
    justifyContent: 'flex-start',
    alignItems: 'center',
    flexDirection: 'column-reverse',
    transform: 'translateX(-50%)'
  },
  vertical: {
    flexDirection: 'row',
    justifyContent: 'flex-end',
    transform: 'translateY(50%)'
  },
  label: {
    lineHeight: 1,
    flexShrink: 0,
    flexGrow: 0
  },
  tick: {
    flexShrink: 0,
    flexGrow: 0
  }
}))

const PlotLabel = React.memo(({
  label,
  size,
  labelPadding,
  tickLength,
  orientation,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = usePlotLabelStyles(classes)
  return <div className={clsx(className, styles[orientation], styles.root)} data-testid={testID}>
    <Typography
      className={styles.label}
      style={{fontSize: size}}
      noWrap
    >{label}
    </Typography>
    <div style={{flex: `0 0 ${labelPadding}px`}} />
    <PlotTick
      className={styles.label}
      length={tickLength}
      orientation={orientation}
    />
  </div>
})

PlotLabel.propTypes = {
  label: PropTypes.string,
  size: PropTypes.number,
  labelPadding: PropTypes.number,
  tickLength: PropTypes.number,
  orientation: PropTypes.oneOf(['horizontal', 'vertical']),
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

PlotLabel.defaultProps = {
  size: 6,
  labelPadding: 3,
  orientation: 'horizontal'
}

export default PlotLabel
