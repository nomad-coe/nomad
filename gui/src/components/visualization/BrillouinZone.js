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
import React, { useState, useEffect, useRef, useCallback } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Typography } from '@material-ui/core'
import {
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay
} from '@material-ui/icons'
import { BrillouinZoneViewer } from '@lauri-codes/materia'
import Floatable from './Floatable'
import NoData from './NoData'
import Placeholder from '../visualization/Placeholder'
import { scale, distance } from '../../utils'
import { withErrorHandler, withWebGLErrorHandler } from '../ErrorHandler'
import { Actions, Action } from '../Actions'

/**
 * Interactive 3D Brillouin zone viewer based on the 'materia'-library.
 */
export const bzError = 'Could not load Brillouin zone.'
const fitMargin = 0.06
const useStyles = makeStyles((theme) => {
  return {
    root: {},
    column: {
      display: 'flex',
      width: '100%',
      height: '100%',
      flexDirection: 'column'
    },
    header: {
      paddingRight: theme.spacing(1),
      display: 'flex',
      flexDirection: 'row',
      zIndex: 1
    },
    canvas: {
      flex: 1,
      zIndex: 0,
      minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
      marginBottom: theme.spacing(2)
    }
  }
})
const BrillouinZone = React.memo(({
  className,
  classes,
  data,
  captureName,
  'data-testid': testID
}) => {
  const [fullscreen, setFullscreen] = useState(false)
  const [loading, setLoading] = useState(true)
  const refViewer = useRef(null)
  const styles = useStyles(classes)

  // In order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const refCanvas = useCallback(node => {
    if (node === null) {
      return
    }
    if (refViewer.current === null) {
      return
    }
    refViewer.current.changeHostElement(node)
    refViewer.current.fit(fitMargin)
    refViewer.current.render()
  }, [])

  // Run only on first render to initialize the viewer.
  const theme = useTheme()
  useEffect(() => {
    if (!data || refViewer.current) {
      return
    }
    refViewer.current = new BrillouinZoneViewer(undefined)
  }, [data])

  // Called only on first render to load the given structure.
  useEffect(() => {
    if (!data) {
      return
    }

    // Format data from section_k_band into the form used by the viewer
    const basis = scale(data.reciprocal_cell, 1E-10)
    const kpoints = []
    const segments = []
    const finalData = {
      basis: basis,
      kpoints: kpoints,
      segments: segments
    }
    let previousPoint
    let segment = []
    for (const seg of data.segment) {
      const labels = seg.endpoints_labels
      const start = seg.kpoints[0]
      const end = seg.kpoints.slice(-1)[0]
      if (!previousPoint || (previousPoint && distance(start, previousPoint) >= 1e-8)) {
        // Push old segment
        if (segment.length > 0) {
          segments.push(segment)
        }
        // Start new segment
        segment = []
        segment.push(start)
        if (labels !== undefined) {
          kpoints.push([labels[0], start])
        }
      }
      segment.push(end)
      if (labels !== undefined) {
        kpoints.push([labels[1], end])
      }
      previousPoint = end
    }
    // Push last segment
    segments.push(segment)

    const options = {
      basis: {
        font: 'Titillium Web,sans-serif',
        offset: 0.025,
        size: 0.04,
        a: { color: '#f44336' },
        b: { color: '#4caf50' },
        c: { color: '#5c6bc0' }
      },
      segments: {
        color: theme.palette.primary.main
      },
      faces: {
        outline: {
          width: 0.002
        }
      },
      kpoints: {
        label: {
          color: theme.palette.primary.main,
          font: 'Titillium Web,sans-serif',
          size: 0.035
        },
        point: {
          color: theme.palette.primary.main,
          size: 0.01
        }
      }
    }

    // Load data into viewer
    refViewer.current.load(finalData, options)
    refViewer.current.setRotation([1, 0, 0, 0])
    refViewer.current.align([
      ['up', 'a'],
      ['front', 'segments']
    ])
    refViewer.current.rotate([
      [0, 1, 0, 45],
      [1, 0, 0, 25]
    ])
    refViewer.current.controls({resetOnDoubleClick: false, pan: {enabled: false}})
    refViewer.current.fit(fitMargin)
    refViewer.current.render()
    refViewer.current.saveCameraReset()
    setLoading(false)
  }, [data, theme])

  const toggleFullscreen = useCallback(() => {
    setFullscreen(!fullscreen)
  }, [fullscreen])

  const takeScreencapture = useCallback(() => {
    refViewer.current.takeScreenShot(captureName)
  }, [captureName])

  const handleReset = useCallback(() => {
    refViewer.current.resetCamera()
    refViewer.current.fit(fitMargin)
    refViewer.current.render()
  }, [])

  // If data is set explicitly to false, we show the NoData component.
  if (data === false) {
    return <NoData className={clsx(className, styles.root)}/>
  }

  if (loading) {
    return <Placeholder
      variant="rect"
      className={clsx(className, styles.root)}
      data-testid={`${testID}-placeholder`}
    />
  }

  return <Floatable
    data-testid={testID}
    className={clsx(styles.root, className)}
    float={fullscreen}
    onFloat={toggleFullscreen}
  >
    <div className={styles.column}>
      {fullscreen && <Typography variant="h6">Brillouin zone</Typography>}
      <div className={styles.canvas} ref={refCanvas}/>
      <div className={styles.header}>
        <Actions>
          <Action tooltip='Reset view' onClick={handleReset}>
            <Replay/>
          </Action>
          <Action tooltip='Toggle fullscreen' onClick={toggleFullscreen}>
            {fullscreen ? <FullscreenExit/> : <Fullscreen/>}
          </Action>
          <Action tooltip='Capture image' onClick={takeScreencapture}>
            <CameraAlt/>
          </Action>
        </Actions>
      </div>
    </div>
  </Floatable>
})

BrillouinZone.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.shape({
      reciprocal_cell: PropTypes.array.isRequired, // Reciprocal cell in SI units
      segment: PropTypes.array.isRequired // Array of section_k_band_segments in SI units
    }),
    PropTypes.oneOf([false, undefined]) // False for NoData, undefined for Placeholder
  ]),
  captureName: PropTypes.string, // Name of the file that the user can download
  classes: PropTypes.object,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}
BrillouinZone.defaultProps = {
  captureName: 'brillouin_zone'
}

export default withWebGLErrorHandler(withErrorHandler(bzError)(BrillouinZone))
