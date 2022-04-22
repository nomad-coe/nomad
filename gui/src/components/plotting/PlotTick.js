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
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'

/**
 * Tick used in plots.
 */
const usePlotTickStyles = makeStyles(theme => ({
  root: {
    backgroundColor: theme.palette.grey[500]
  }
}))
const PlotTick = React.memo(({
  length,
  thickness,
  orientation,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = usePlotTickStyles(classes)
  return <div
    className={clsx(className, styles.root)}
    style={orientation === 'vertical'
      ? {width: length, height: thickness}
      : {width: thickness, height: length}
    }
    data-testid={testID}
  />
})

PlotTick.propTypes = {
  length: PropTypes.number,
  thickness: PropTypes.number,
  orientation: PropTypes.oneOf(['horizontal', 'vertical']),
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

PlotTick.defaultProps = {
  length: 6,
  thickness: 1
}

export default PlotTick
