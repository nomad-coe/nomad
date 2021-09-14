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
import { Typography } from '@material-ui/core/'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'

/**
 * A rectangular bar displaying the relative occurence of a specific value. Uses
 * an animated CSS transform to react to changes in the bar size. Notice that
 * using a transform is much more performant than animating the width.
 */
const useStyles = makeStyles(theme => ({
  root: {
    zIndex: -100
  },
  container: {
    position: 'relative',
    width: '100%',
    height: '100%'
  },
  rectangle: {
    backgroundColor: theme.palette.secondary.light,
    width: '100%',
    height: '100%',
    '-webkit-transform': 'none',
    transform: 'none',
    transition: 'transform 250ms',
    transformOrigin: 'center left',
    willChange: 'transform'
  },
  value: {
    position: 'absolute',
    right: theme.spacing(1),
    top: 0,
    bottom: 0
  }
}))
const InputStatisticsBar = React.memo(({
  max,
  value,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = useStyles()
  const scale = (value && max) ? value / max : 0

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div className={styles.container}>
      <div className={styles.rectangle} style={{transform: `scaleX(${scale})`}}></div>
      <Typography className={styles.value}>{value || ''}</Typography>
    </div>
  </div>
})

InputStatisticsBar.propTypes = {
  max: PropTypes.number,
  value: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default InputStatisticsBar
