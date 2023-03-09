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
import React, { useCallback, useMemo } from 'react'
import clsx from 'clsx'
import { clamp } from 'lodash'
import { Tooltip } from '@material-ui/core/'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'

/**
 * Bar within a plot.
 */
const usePlotAxisStyles = makeStyles(theme => ({
  root: {
    position: 'absolute',
    cursor: 'pointer',
    height: '100%',
    paddingLeft: '1px',
    paddingRight: '1px'
  },
  container: {
    position: 'relative',
    width: '100%',
    height: '100%'
  },
  rectangle: {
    width: '100%',
    height: '100%',
    '-webkit-transform': 'none',
    transform: 'none',
    transition: 'transform 250ms',
    transformOrigin: 'bottom left',
    willChange: 'transform',
    backgroundColor: theme.palette.primary.veryLight,
    position: 'relatives'
  },
  highlight: {
    position: 'absolute',
    top: 0,
    bottom: 0,
    backgroundColor: theme.palette.secondary.veryLight
  }
}))

const PlotBar = React.memo(({
  startX,
  endX,
  startY,
  endY,
  selected,
  tooltip,
  onClick,
  className,
  classes,
  'data-testid': testID}) => {
  const styles = usePlotAxisStyles(classes)

  const handleClick = useCallback((event) => {
    onClick && onClick(event)
  }, [onClick])

  const highlightStyle = useMemo(() => {
    if (selected === false) {
      return {visibility: 'hidden'}
    }
    if (selected === true) {
      return {
        left: '0%',
        right: '0%'
      }
    }
    return {
      left: `${clamp(selected[0] * 100, 0, 100)}%`,
      right: `${clamp((1 - selected[1]) * 100, 0, 100)}%`
    }
  }, [selected])

  return <div
    onClick={handleClick}
    className={clsx(className, styles.root)}
    style={{left: `${startX * 100}%`, right: `${(1 - endX) * 100}%`}}
    data-testid={testID}
  >
    <Tooltip placement="bottom" enterDelay={0} leaveDelay={0} title={tooltip || ''} arrow>
      <div className={styles.container}>
        <div className={styles.rectangle} style={{transform: `scaleY(${endY - startY})`}}>
          <div
            className={styles.highlight}
            style={highlightStyle}/>
        </div>
      </div>
    </Tooltip>
  </div>
})

PlotBar.propTypes = {
  startX: PropTypes.number,
  endX: PropTypes.number,
  startY: PropTypes.number,
  endY: PropTypes.number,
  /* The selection as a boolean, or as a range. The range is given as two
   * percentages indicate the start and end of the selection. */
  selected: PropTypes.oneOfType([PropTypes.bool, PropTypes.arrayOf(PropTypes.number)]),
  tooltip: PropTypes.string,
  onClick: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

PlotBar.defaultProps = {
}

export default PlotBar
