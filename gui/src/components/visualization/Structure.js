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
import { makeStyles } from '@material-ui/core/styles'
import {
  Box,
  Checkbox,
  Menu,
  MenuItem,
  IconButton,
  Tooltip,
  Typography,
  FormControlLabel
} from '@material-ui/core'
import {
  MoreVert,
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay
} from '@material-ui/icons'
import { StructureViewer } from '@lauri-codes/materia'
import Floatable from './Floatable'

export default function Structure({className, classes, system, options, viewer, captureName, aspectRatio, positionsOnly}) {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showBonds, setShowBonds] = useState(true)
  const [showLatticeConstants, setShowLatticeConstants] = useState(true)
  const [showCell, setShowCell] = useState(true)

  // Variables
  const open = Boolean(anchorEl)
  const refViewer = useRef(null)
  const refCanvas = useRef(null)

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
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
      iconButton: {
        backgroundColor: 'white'
      }
    }
  })
  const style = useStyles(classes)

  // In order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const measuredRef = useCallback(node => {
    refCanvas.current = node
    if (node === null) {
      return
    }
    if (refViewer.current === null) {
      return
    }
    refViewer.current.changeHostElement(node, true, true)
  }, [])

  // Run only on first render to initialize the viewer. See the viewer
  // documentation for details on the meaning of different options:
  // https://lauri-codes.github.io/materia/viewers/structureviewer
  useEffect(() => {
    let viewerOptions
    if (options === undefined) {
      viewerOptions = {
        view: {
          autoResize: true,
          autoFit: true,
          fitMargin: 0.5
        },
        bonds: {
          enabled: true
        },
        latticeConstants: {
          size: 0.7,
          font: 'Titillium Web,sans-serif',
          a: {color: '#f44336'},
          b: {color: '#4caf50'},
          c: {color: '#5c6bc0'}
        },
        controls: {
          enableZoom: true,
          enablePan: true,
          enableRotate: true
        },
        renderer: {
          backgroundColor: ['#ffffff', 1],
          shadows: {
            enabled: false
          }
        }
      }
    } else {
      viewerOptions = options
    }
    if (viewer === undefined) {
      refViewer.current = new StructureViewer(undefined, viewerOptions)
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
    if (system === undefined) {
      return
    }

    if (positionsOnly) {
      refViewer.current.setPositions(system.positions)
      return
    }

    // Systems with cell are centered on the cell center and orientation is defined
    // by the cell vectors.
    let cell = system.cell
    if (cell !== undefined) {
      refViewer.current.setOptions({layout: {
        viewCenter: 'COC',
        viewRotation: {
          alignments: [
            ['up', 'c'],
            ['right', 'b']
          ],
          rotations: [
            [0, 1, 0, 60],
            [1, 0, 0, 30]
          ]
        }
      }}, false, false)
    // Systems without cell are centered on the center of positions
    } else {
      refViewer.current.setOptions({layout: {
        viewCenter: 'COP',
        viewRotation: {
          rotations: [
            [0, 1, 0, 60],
            [1, 0, 0, 30]
          ]
        }
      }}, false, false)
    }
    refViewer.current.load(system)
    refViewer.current.fitToCanvas()
    refViewer.current.saveReset()
    refViewer.current.reset()
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Viewer settings
  useEffect(() => {
    refViewer.current.setOptions({bonds: {enabled: showBonds}})
  }, [showBonds])

  useEffect(() => {
    refViewer.current.setOptions({latticeConstants: {enabled: showLatticeConstants}})
  }, [showLatticeConstants])

  useEffect(() => {
    refViewer.current.setOptions({cell: {enabled: showCell}})
  }, [showCell])

  // Memoized callbacks
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])

  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

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

  const content = <Box className={style.container}>
    <div className={style.header}>
      {fullscreen && <Typography variant="h6">Structure</Typography>}
      <div className={style.spacer}></div>
      <Tooltip title="Reset view">
        <IconButton className={style.iconButton} onClick={handleReset}>
          <Replay />
        </IconButton>
      </Tooltip>
      <Tooltip
        title="Toggle fullscreen">
        <IconButton className={style.iconButton} onClick={toggleFullscreen}>
          {fullscreen ? <FullscreenExit /> : <Fullscreen />}
        </IconButton>
      </Tooltip>
      <Tooltip title="Capture image">
        <IconButton className={style.iconButton} onClick={takeScreencapture}>
          <CameraAlt />
        </IconButton>
      </Tooltip>
      <Tooltip title="Options">
        <IconButton className={style.iconButton} onClick={openMenu}>
          <MoreVert />
        </IconButton>
      </Tooltip>
      <Menu
        id='settings-menu'
        anchorEl={anchorEl}
        getContentAnchorEl={null}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
        transformOrigin={{ vertical: 'top', horizontal: 'right' }}
        keepMounted
        open={open}
        onClose={closeMenu}
      >
        <MenuItem key='show-bonds'>
          <FormControlLabel
            control={
              <Checkbox
                checked={showBonds}
                onChange={(event) => { setShowBonds(!showBonds) }}
                color='primary'
              />
            }
            label='Show bonds'
          />
        </MenuItem>
        <MenuItem key='show-axis'>
          <FormControlLabel
            control={
              <Checkbox
                checked={showLatticeConstants}
                onChange={(event) => { setShowLatticeConstants(!showLatticeConstants) }}
                color='primary'
              />
            }
            label='Show lattice constants'
          />
        </MenuItem>
        <MenuItem key='show-cell'>
          <FormControlLabel
            control={
              <Checkbox
                checked={showCell}
                onChange={(event) => { setShowCell(!showCell) }}
                color='primary'
              />
            }
            label='Show simulation cell'
          />
        </MenuItem>
      </Menu>
    </div>
    <div className={style.viewerCanvas} ref={measuredRef}></div>
  </Box>

  return (
    <Floatable float={fullscreen} onFloat={toggleFullscreen} aspectRatio={aspectRatio}>
      {content}
    </Floatable>
  )
}

Structure.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  viewer: PropTypes.object, // Optional shared viewer instance.
  system: PropTypes.object, // The system to display as section_system
  options: PropTypes.object, // Viewer options
  captureName: PropTypes.string, // Name of the file that the user can download
  aspectRatio: PropTypes.number, // Fixed aspect ratio for the viewer canvas
  positionsOnly: PropTypes.bool // Whether to update only positions. This is much faster than loading the entire structure.
}
Structure.defaultProps = {
  aspectRatio: 4 / 3,
  captureName: 'structure'
}
