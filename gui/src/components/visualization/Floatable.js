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
import React, { useEffect, useCallback, useMemo, useRef } from 'react'
import { useResizeDetector } from 'react-resize-detector'
import { makeStyles } from '@material-ui/core/styles'
import {
  Button,
  Dialog,
  DialogContent,
  DialogActions,
  Box,
  Paper,
  Backdrop
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import { useWindowSize } from '../../hooks'

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

/**
 * Component that wraps it's children in a container that can be 'floated', i.e.
 * is positioned relative to the viewport and is above all other elements.
 */
const actionsHeight = 52.5 // Size of the actions element that is shown when floating
const useNoReparentStyles = makeStyles((theme) => {
  return {
    root: {
      top: '50%',
      left: '50%',
      margin: 0,
      padding: 0,
      maxWidth: 'none',
      maxHeight: 'none',
      boxSizing: 'border-box'
    },
    paper: {
      display: 'flex',
      flexDirection: 'column',
      width: '100%',
      height: '100%',
      boxSizing: 'border-box'
    },
    full: {
      width: '100%',
      height: '100%',
      boxSizing: 'border-box'
    },
    dialogActions: {
      height: actionsHeight,
      boxSizing: 'border-box'
    }
  }
})
export const FloatableNoReparent = React.memo(({
  className,
  float,
  children,
  onFloat,
  onResize,
  maxWidth,
  margin
}) => {
  // The ratio is calculated dynamically from the final size
  const styles = useNoReparentStyles()
  const { height, width, ref } = useResizeDetector()
  const dim = useRef({width: 800, height: 600})
  const windowSize = useWindowSize()

  // This fixes the aspect ratio of the component based on the non-floating
  // size.
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

  const handleFloat = useCallback(() => {
    onFloat && onFloat()
  }, [onFloat])

  // Calculate the final floated size
  const size = useMemo(() => {
    // Calculate the fullscreen size
    const marginAbsolute = margin * (windowSize.windowHeight)
    const windowHeight = windowSize.windowHeight - marginAbsolute
    const windowWidth = Math.min(windowSize.windowWidth, maxWidth) - marginAbsolute
    const windowRatio = windowWidth / windowHeight
    const ratio = dim.current.width / dim.current.height
    let width
    let height
    if (windowRatio > ratio) {
      width = (windowHeight) * ratio
      height = windowHeight + actionsHeight
    } else {
      width = windowWidth
      height = (windowWidth) / ratio + actionsHeight
    }
    return {width, height}
  }, [windowSize, maxWidth, margin])

  return <div
    ref={ref}
    className={clsx(styles.root, className)}
    style={{
      position: float ? 'fixed' : 'static',
      zIndex: float ? 10000 : 1,
      width: float ? size.width : '100%',
      height: float ? size.height : '100%',
      marginTop: float ? -0.5 * size.height : 'unset',
      marginLeft: float ? -0.5 * size.width : 'unset'
    }}
  >
    <Paper elevation={float ? 1 : 0} className={styles.paper}>
      <Box padding={float ? 3 : 0} paddingBottom={float ? 1 : 0} className={styles.full}>
        <div></div>
        {children}
      </Box>
      {float &&
        <DialogActions className={styles.dialogActions}>
          <Button onClick={handleFloat}>
            Close
          </Button>
        </DialogActions>
      }
    </Paper>
    <Backdrop className={styles.backdrop} open={float} onClick={handleFloat} />
  </div>
})

FloatableNoReparent.propTypes = {
  float: PropTypes.bool,
  onFloat: PropTypes.func,
  onResize: PropTypes.func,
  children: PropTypes.node,
  className: PropTypes.string,
  maxWidth: PropTypes.number,
  margin: PropTypes.number
}

FloatableNoReparent.defaultProps = {
  float: false,
  maxWidth: 1920,
  margin: 0.1
}
