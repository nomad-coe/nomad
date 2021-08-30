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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Box,
  Typography
} from '@material-ui/core'
import {
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay
} from '@material-ui/icons'
import { BrillouinZoneViewer } from '@lauri-codes/materia'
import Floatable from './Floatable'
import Placeholder from '../visualization/Placeholder'
import { scale, distance } from '../../utils'
import { withErrorHandler, withWebGLErrorHandler } from '../ErrorHandler'
import { Actions, Action } from '../Actions'
import clsx from 'clsx'

/**
 * Interactive 3D Brillouin zone viewer based on the 'materia'-library.
 */
const BrillouinZone = React.memo(({
  className,
  classes,
  options,
  viewer,
  data,
  captureName,
  aspectRatio,
  'data-testid': testID
}) => {
  // States
  const [fullscreen, setFullscreen] = useState(false)
  const [loading, setLoading] = useState(true)

  // Refs
  const refViewer = useRef(null)

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      root: {
      },
      container: {
        display: 'flex',
        width: '100%',
        height: '100%',
        flexDirection: 'column',
        backgroundColor: 'white'
      },
      header: {
        paddingRight: theme.spacing(1),
        display: 'flex',
        flexDirection: 'row',
        zIndex: 1
      },
      spacer: {
        flex: 1
      },
      viewerCanvas: {
        flex: 1,
        zIndex: 0,
        minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
        marginBottom: theme.spacing(2)
      },
      errorMessage: {
        flex: '0 0 70%',
        color: '#aaa',
        textAlign: 'center'
      },
      iconButton: {
        backgroundColor: 'white'
      }
    }
  })
  let style = useStyles(classes)

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
    refViewer.current.changeHostElement(node, true, true)
    refCanvas.current = node
  }, [])

  // Run only on first render to initialize the viewer. See the viewer
  // documentation for details on the meaning of different options:
  // https://nomad-coe.github.io/materia/viewers/brillouinzoneviewer
  const theme = useTheme()
  useEffect(() => {
    let viewerOptions
    if (options === undefined) {
      viewerOptions = {
        view: {
          autoResize: false,
          autoFit: true,
          fitMargin: 0.06
        },
        layout: {
          viewRotation: {
            alignments: [
              ['up', 'a'],
              ['front', 'segments']
            ],
            rotations: [
              [0, 1, 0, 45],
              [1, 0, 0, 25]
            ]
          }
        },
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
    } else {
      viewerOptions = options
    }

    if (viewer === undefined) {
      refViewer.current = new BrillouinZoneViewer(undefined, viewerOptions)
    } else {
      refViewer.current = viewer
      refViewer.current.setOptions(viewerOptions, false, false)
    }
    if (refCanvas.current !== null) {
      refViewer.current.changeHostElement(refCanvas.current, false, false)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

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
    for (let seg of data.segments) {
      let labels = [seg.kpoints_labels[0], seg.kpoints_labels.pop()]
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

    // Load data into viewer
    refViewer.current.load(finalData)
    refViewer.current.fitToCanvas()
    refViewer.current.saveReset()
    refViewer.current.reset()
    setLoading(false)
  }, [data])

  const toggleFullscreen = useCallback(() => {
    setFullscreen(!fullscreen)
  }, [fullscreen])

  const takeScreencapture = useCallback(() => {
    refViewer.current.takeScreenShot(captureName)
  }, [captureName])

  const handleReset = useCallback(() => {
    refViewer.current.reset()
    refViewer.current.fitToCanvas()
    refViewer.current.render()
  }, [])

  if (loading) {
    return <Placeholder
      variant="rect"
      aspectRatio={aspectRatio}
      data-testid={`${testID}-placeholder`}
    />
  }

  const content = <Box className={style.container}>
    {fullscreen && <Typography variant="h6">Brillouin zone</Typography>}
    <div className={style.viewerCanvas} ref={refCanvas}></div>
    <div className={style.header}>
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
  </Box>

  return <Box className={clsx(style.root, className)} >
    <Floatable float={fullscreen} onFloat={toggleFullscreen} aspectRatio={aspectRatio}>
      {content}
    </Floatable>
  </Box>
})

BrillouinZone.propTypes = {
  viewer: PropTypes.object, // Optional shared viewer instance.
  data: PropTypes.shape({
    reciprocal_cell: PropTypes.array.isRequired, // Reciprocal cell in SI units
    segments: PropTypes.array.isRequired // Array of section_k_band_segments in SI units
  }),
  options: PropTypes.object, // Viewer options
  captureName: PropTypes.string, // Name of the file that the user can download
  aspectRatio: PropTypes.number, // Fixed aspect ratio for the viewer canvas
  classes: PropTypes.object,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}
BrillouinZone.defaultProps = {
  aspectRatio: 4 / 3,
  captureName: 'brillouin_zone'
}

export default withWebGLErrorHandler(withErrorHandler(BrillouinZone, 'Could not load Brillouin zone.'))
