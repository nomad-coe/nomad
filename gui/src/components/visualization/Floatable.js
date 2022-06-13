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
import React, { useEffect, useRef } from 'react'
import { useResizeDetector } from 'react-resize-detector'
import { makeStyles } from '@material-ui/core/styles'
import {
  Button,
  Dialog,
  DialogContent,
  DialogActions
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'

/**
 * Component that wraps it's children in a container that can be 'floated',
 * i.e. displayed on an html element that is positioned relative to the
 * viewport and is above all other elements.
 */
const useStyles = makeStyles((theme) => {
  const actionsHeight = 52.5 // Size of the actions element that is shown when floating
  const padding = 20 // Padding arount the floating object
  return {
    root: {
      width: '100%',
      height: '100%'
    },
    dialogContent: {
      boxSizing: 'border-box',
      padding: padding
    },
    'dialogContent:first-child': {
      paddingTop: `${padding} !important`
    },
    dialogActions: {
      height: actionsHeight,
      boxSizing: 'border-box'
    }
  }
})
export default function Floatable({
  className,
  classes,
  style,
  float,
  children,
  onFloat,
  onResize
}) {
  // The ratio is calculated dynamically from the final size
  const { height, width, ref } = useResizeDetector()
  const dim = useRef({width: 800, height: 600})

  // When the figure is floated, fix the width/height.
  useEffect(() => {
    if (!float) {
      dim.current = {
        width: width,
        height: height
      }
    }
    if (onResize && !isNil(width) && !isNil(height)) {
      onResize(width, height)
    }
  }, [float, width, height, onResize])

  // Dynamic styles
  const useDynamicStyles = makeStyles((theme) => {
    // Calculate the fullscreen size
    const actionsHeight = 52.5 // Size of the actions element that is shown when floating
    const padding = 20 // Padding arount the floating object
    const maxWidth = 1280 // Maximum width for the floating window
    const margin = 0.1 * (window.innerHeight) + actionsHeight
    const windowHeight = window.innerHeight - margin
    const windowWidth = Math.min(window.innerWidth, maxWidth) - margin
    const windowRatio = windowWidth / windowHeight
    const ratio = dim.current.width / dim.current.height
    let width
    let height
    if (windowRatio > ratio) {
      width = (windowHeight) * ratio + 2 * padding + 'px'
      height = windowHeight + actionsHeight + 2 * padding
    } else {
      width = windowWidth + 2 * padding
      height = (windowWidth) / ratio + actionsHeight + 2 * padding + 'px'
    }

    return {
      dialogRoot: {
        width: width,
        height: height,
        margin: 0,
        padding: 0,
        maxWidth: 'none',
        maxHeight: 'none',
        boxSizing: 'border-box'
      }
    }
  })
  const styles = useStyles({classes: classes})
  const dynamicStyles = useDynamicStyles({classes: classes})

  return <div ref={ref} className={clsx(styles.root, className)} style={style}>
    {float
      ? <div style={{width: dim.current.width, height: dim.current.height}} />
      : children
    }
    <Dialog fullWidth={false} maxWidth={false} open={float}
      classes={{paper: dynamicStyles.dialogRoot}}
    >
      <DialogContent className={[styles.dialogContent, styles['dialogContent:first-child']].join('_')}>
        {float ? children : ''}
      </DialogContent>
      <DialogActions className={styles.dialogActions}>
        <Button onClick={() => onFloat(float)}>
          Close
        </Button>
      </DialogActions>
    </Dialog>
  </div>
}

Floatable.propTypes = {
  /**
   * Whether this component should be in floating mode or not.
   */
  float: PropTypes.bool.isRequired,
  /**
   * Callback that is called whenever this component requests a change in the
   * float property. The callback accepts one parameter: 'float' that is a
   * boolean indicating the current float status.
   */
  onFloat: PropTypes.func,
  onResize: PropTypes.func,
  children: PropTypes.node,
  className: PropTypes.string,
  classes: PropTypes.object,
  style: PropTypes.object
}
Floatable.defaultProps = {
  float: false
}
