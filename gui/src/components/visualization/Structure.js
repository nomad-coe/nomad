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
import React, { useState, useEffect, useRef, useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import { makeStyles, fade } from '@material-ui/core/styles'
import {
  Box,
  Checkbox,
  Menu,
  MenuItem,
  Typography,
  FormControlLabel
} from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup } from '@material-ui/lab'
import {
  MoreVert,
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay
} from '@material-ui/icons'
import { StructureViewer } from '@lauri-codes/materia'
import Floatable from './Floatable'
import Placeholder from '../visualization/Placeholder'
import Actions from '../Actions'
import { mergeObjects } from '../../utils'
import { withErrorHandler, ErrorCard } from '../ErrorHandler'
import _ from 'lodash'
import clsx from 'clsx'

/**
 * Used to show atomistic systems in an interactive 3D viewer based on the
 * 'materia'-library.
 */
export const Structure = withErrorHandler(({
  className,
  classes,
  system,
  systems,
  options,
  viewer,
  captureName,
  aspectRatio,
  positionsOnly,
  sizeLimit,
  positionsSubject}
) => {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showBonds, setShowBonds] = useState(true)
  const [showLatticeConstants, setShowLatticeConstants] = useState(true)
  const [showCell, setShowCell] = useState(true)
  const [wrap, setWrap] = useState(true)
  const [showPrompt, setShowPrompt] = useState(false)
  const [accepted, setAccepted] = useState(false)
  const [nAtoms, setNAtoms] = useState(false)
  const [loading, setLoading] = useState(true)
  const [shownSystem, setShownSystem] = useState(null)
  const [finalSystem, setFinalSystem] = useState(system)

  // Variables
  const open = Boolean(anchorEl)
  const refViewer = useRef(null)

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      root: {},
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
      toggles: {
        marginBottom: theme.spacing(1),
        height: '2rem'
      },
      toggle: {
        color: fade(theme.palette.action.active, 0.87)
      },
      selected: {
        '&$selected': {
          color: fade(theme.palette.action.active, 0.87)
        }
      },
      title: {
        marginBottom: theme.spacing(1)
      },
      viewerCanvas: {
        flex: 1,
        zIndex: 0,
        minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
        marginBottom: theme.spacing(1)
      }
    }
  })
  const styles = useStyles(classes)

  useEffect(() => {
    setFinalSystem(system)
  }, [system])

  useEffect(() => {
    if (systems) {
      const firstSystem = Object.keys(systems)[0]
      setFinalSystem(systems[firstSystem])
      setShownSystem(firstSystem)
    }
  }, [systems])

  // In order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const refCanvas = useCallback(node => {
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
  // https://nomad-coe.github.io/materia/viewers/structureviewer
  useEffect(() => {
    let viewerOptions
    if (options === undefined) {
      viewerOptions = {
        view: {
          autoResize: true,
          autoFit: true,
          fitMargin: 0.6
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
      viewerOptions = mergeObjects(options, viewerOptions)
    }
    if (viewer === undefined) {
      refViewer.current = new StructureViewer(undefined, viewerOptions)
    } else {
      refViewer.current = viewer
      refViewer.current.setOptions(viewerOptions, false, false)
    }
    if (refCanvas.current) {
      refViewer.current.changeHostElement(refCanvas.current, false, false)
    }
    if (positionsSubject && refViewer.current) {
      positionsSubject.subscribe((positions) => {
        refViewer.current.setPositions(positions)
      })
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [positionsSubject])

  const loadSystem = useCallback((system, refViewer) => {
    // If the cell all zeroes, positions are assumed to be cartesian.
    if (system.cell !== undefined) {
      if (_.sum(_.flattenDeep(system.cell)) === 0) {
        system.cell = undefined
      }
    }
    // Systems with non-empty cell are centered on the cell center and
    // orientation is defined by the cell vectors.
    if (system.cell !== undefined) {
      refViewer.current.setOptions({layout: {
        viewCenter: 'COC',
        periodicity: 'wrap',
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
    setLoading(false)
  }, [])

  // Called whenever the given system changes. If positionsOnly is true, only
  // updates the positions. Otherwise reloads the entire structure.
  useEffect(() => {
    if (!finalSystem) {
      return
    }

    if (!accepted) {
      const nAtoms = finalSystem.positions.length
      setNAtoms(nAtoms)
      if (nAtoms > 300) {
        setShowPrompt(true)
        return
      }
    }

    if (positionsOnly && !!(refViewer?.current?.structure)) {
      refViewer.current.setPositions(finalSystem.positions)
      setLoading(false)
      return
    }

    loadSystem(finalSystem, refViewer)
  }, [finalSystem, positionsOnly, loadSystem, accepted])

  // Viewer settings
  useEffect(() => {
    refViewer.current.setOptions({bonds: {enabled: showBonds}})
  }, [showBonds])

  useEffect(() => {
    refViewer.current.setOptions({latticeConstants: {enabled: showLatticeConstants}})
  }, [showLatticeConstants])

  useEffect(() => {
    refViewer.current.setOptions({layout: {periodicity: wrap ? 'wrap' : 'none'}})
  }, [wrap])

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

  const structureToggles = useMemo(() => {
    if (systems) {
      const toggles = []
      for (let key in systems) {
        toggles.push(<ToggleButton
          key={key}
          value={key}
          aria-label={key}
          classes={{root: styles.toggle, selected: styles.selected}}
        >{key}</ToggleButton>)
      }
      return toggles
    }
    return null
  }, [systems, styles])

  // Enforce at least one structure view option
  const handleStructureChange = (event, value) => {
    if (value !== null) {
      setShownSystem(value)
      setFinalSystem(systems[value])
    }
  }

  if (showPrompt) {
    return <ErrorCard
      message={`Visualization is by default disabled for systems with more than ${sizeLimit} atoms. Do you wish to enable visualization for this system with ${nAtoms} atoms?`}
      className={styles.error}
      actions={[{label: 'Yes', onClick: e => { setShowPrompt(false); loadSystem(finalSystem, refViewer); setAccepted(true) }}]}
    >
    </ErrorCard>
  }
  if (loading) {
    return <Placeholder variant="rect" aspectRatio={aspectRatio}></Placeholder>
  }

  // List of actionable buttons for the viewer
  const actions = [
    {tooltip: 'Reset view', onClick: handleReset, content: <Replay/>},
    {tooltip: 'Toggle fullscreen', onClick: toggleFullscreen, content: fullscreen ? <FullscreenExit/> : <Fullscreen/>},
    {tooltip: 'Capture image', onClick: takeScreencapture, content: <CameraAlt/>},
    {tooltip: 'Options', onClick: openMenu, content: <MoreVert/>}
  ]

  const content = <Box className={styles.container}>
    {fullscreen && <Typography className={styles.title} variant="h6">Structure</Typography>}
    {systems && <ToggleButtonGroup
      className={styles.toggles}
      size="small"
      exclusive
      value={shownSystem}
      onChange={handleStructureChange}
    >
      {structureToggles}
    </ToggleButtonGroup>
    }
    <div className={styles.viewerCanvas} ref={refCanvas}></div>
    <div className={styles.header}>
      <Actions actions={actions}></Actions>
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
        <MenuItem key='wrap'>
          <FormControlLabel
            control={
              <Checkbox
                checked={wrap}
                onChange={(event) => { setWrap(!wrap) }}
                color='primary'
              />
            }
            label='Wrap positions'
          />
        </MenuItem>
      </Menu>
    </div>
  </Box>
  return <Box className={clsx(styles.root, className)} >
    <Floatable float={fullscreen} onFloat={toggleFullscreen} aspectRatio={aspectRatio}>
      {content}
    </Floatable>
  </Box>
}, 'Could not load structure.')

Structure.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  viewer: PropTypes.object, // Optional shared viewer instance.
  system: PropTypes.object, // The system to display in the native materia-format
  systems: PropTypes.object, // Set of systems that can be switched
  options: PropTypes.object, // Viewer options
  captureName: PropTypes.string, // Name of the file that the user can download
  aspectRatio: PropTypes.number, // Fixed aspect ratio for the viewer canvas
  positionsOnly: PropTypes.bool, // Whether to update only positions. This is much faster than loading the entire structure.
  sizeLimit: PropTypes.number, // Maximum system size before a prompt is shown
  /**
   * A RxJS Subject for efficient, non-persistent, position changes that bypass
   * rendering of the component. Should send messages that contain the new
   * atomic positions as a list.
  */
  positionsSubject: PropTypes.any
}
Structure.defaultProps = {
  aspectRatio: 4 / 3,
  captureName: 'structure',
  sizeLimit: 300
}
