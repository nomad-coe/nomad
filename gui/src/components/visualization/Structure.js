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
import { makeStyles, alpha } from '@material-ui/core/styles'
import {
  Checkbox,
  Button,
  Menu,
  MenuItem,
  Typography,
  FormControlLabel
} from '@material-ui/core'
import { Alert } from '@material-ui/lab'
import {
  MoreVert,
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay,
  ViewList
} from '@material-ui/icons'
import { StructureViewer } from '@lauri-codes/materia'
import Floatable from './Floatable'
import NoData from './NoData'
import Placeholder from '../visualization/Placeholder'
import { Actions, Action } from '../Actions'
import { mergeObjects, delay } from '../../utils'
import { withErrorHandler, withWebGLErrorHandler } from '../ErrorHandler'
import { useHistory } from 'react-router-dom'
import _ from 'lodash'
import clsx from 'clsx'

/**
 * Used to show atomistic systems in an interactive 3D viewer based on the
 * 'materia'-library.
 */
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
      display: 'flex',
      flexDirection: 'row',
      zIndex: 1
    },
    toggles: {
      marginBottom: theme.spacing(1),
      height: '2rem'
    },
    toggle: {
      color: alpha(theme.palette.action.active, 0.87)
    },
    selected: {
      '&$selected': {
        color: alpha(theme.palette.action.active, 0.87)
      }
    },
    title: {
      marginBottom: theme.spacing(1)
    },
    canvas: {
      flex: 1,
      zIndex: 0,
      minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
      marginBottom: theme.spacing(1)
    },
    switchContainer: {
      flex: 1,
      zIndex: 0,
      minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
      marginBottom: theme.spacing(1)
    },
    switchPlaceholder: {
      left: theme.spacing(0),
      right: theme.spacing(0),
      top: theme.spacing(0.5),
      bottom: theme.spacing(0.5)
    }
  }
})
const Structure = React.memo(({
  className,
  classes,
  data,
  options,
  materialType,
  structureType,
  m_path,
  viewer,
  captureName,
  positionsOnly,
  sizeLimit,
  bondLimit,
  positionsSubject,
  'data-testid': testID,
  placeHolderStyle,
  noDataStyle}
) => {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showLatticeConstants, setShowLatticeConstants] = useState(true)
  const [showCell, setShowCell] = useState(true)
  const [wrap, setWrap] = useState(materialType === 'bulk')
  const [showPrompt, setShowPrompt] = useState(false)
  const [accepted, setAccepted] = useState(false)
  const [nAtoms, setNAtoms] = useState()
  const [showBonds, setShowBonds] = useState(false)
  const [loading, setLoading] = useState(true)

  // Variables
  const history = useHistory()
  const open = Boolean(anchorEl)
  const refViewer = useRef(null)
  const styles = useStyles(classes)

  useEffect(() => {
    if (data) {
      setLoading(true)
    }
  }, [data])

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

  // Merge custom options with default options. These options can only set once.
  const finalOptions = useMemo(() => {
    const defaultOptions = {
      view: {
        autoResize: false,
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
      layout: {
        periodicity: wrap ? 'wrap' : 'none'
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
    let viewerOptions
    if (options === undefined) {
      viewerOptions = defaultOptions
    } else {
      viewerOptions = mergeObjects(options, defaultOptions)
    }
    return viewerOptions
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Run only on first render to initialize the viewer. See the viewer
  // documentation for details on the meaning of different options:
  // https://nomad-coe.github.io/materia/viewers/structureviewer
  useEffect(() => {
    if (viewer === undefined) {
      refViewer.current = new StructureViewer(undefined, finalOptions)
    } else {
      refViewer.current = viewer
      refViewer.current.setOptions(finalOptions, false, false)
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
    // This function calls is delayed in order to not block the first render of
    // the component. delay hoists the heavy function call outside so that the
    // react rendering can finish before it.
    delay(() => {
      // If the cell is all zeroes, positions are assumed to be cartesian.
      if (system.cell !== undefined) {
        if (_.flattenDeep(system.cell).every(v => v === 0)) {
          system.cell = undefined
        }
      }
      // Determine the orientation and view centering based on material type and
      // the structure type.
      let layout = {
        viewCenter: 'COP',
        viewRotation: {
          alignments: system.cell === undefined
            ? undefined
            : [
              ['up', 'c'],
              ['right', 'b']
            ],
          rotations: [
            [0, 1, 0, 60],
            [1, 0, 0, 30]
          ]
        }
      }
      if (structureType === 'conventional' || structureType === 'primitive') {
        if (materialType === 'bulk') {
          layout.viewCenter = 'COC'
          layout.viewRotation.alignments = [
            ['up', 'c'],
            ['right', 'b']
          ]
        } else if (materialType === '2D') {
          layout = {
            viewCenter: 'COC',
            viewRotation: {
              alignments: [
                ['right', 'a'],
                ['up', 'b']
              ],
              rotations: [
                [1, 0, 0, -60]
              ]
            }
          }
        } else if (materialType === '1D') {
          layout = {
            viewCenter: 'COC',
            viewRotation: {
              alignments: [
                ['right', 'a']
              ],
              rotations: [
                [1, 0, 0, -60],
                [0, 1, 0, 30],
                [1, 0, 0, 30]
              ]
            }
          }
        }
      }
      const opts = {layout}
      refViewer.current.setOptions(opts, false, false)
      refViewer.current.load(system)
      refViewer.current.fitToCanvas()
      refViewer.current.saveReset()
      refViewer.current.reset()
      setLoading(false)
    })
  }, [materialType, structureType])

  // Called whenever the given system changes. If positionsOnly is true, only
  // updates the positions. Otherwise reloads the entire structure.
  useEffect(() => {
    if (!data) {
      return
    }

    if (!accepted) {
      const nAtoms = data.positions.length
      setNAtoms(nAtoms)
      if (nAtoms <= bondLimit) {
        setShowBonds(true)
      }
      if (nAtoms > sizeLimit) {
        setShowPrompt(true)
        return
      }
    }

    if (positionsOnly && !!(refViewer?.current?.structure)) {
      refViewer.current.setPositions(data.positions)
      setLoading(false)
      return
    }

    loadSystem(data, refViewer)
  }, [data, positionsOnly, accepted, loadSystem, sizeLimit, bondLimit])

  // Viewer settings
  useEffect(() => {
    refViewer.current.setOptions({bonds: {enabled: showBonds}})
  }, [showBonds])

  useEffect(() => {
    refViewer.current.setOptions({latticeConstants: {enabled: showLatticeConstants}})
  }, [showLatticeConstants])

  useEffect(() => {
    if (data?.cell) {
      refViewer.current.setOptions({layout: {periodicity: wrap ? 'wrap' : 'none'}})
    }
  }, [wrap, data])

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

  // If data is set explicitly to false, we show the NoData component.
  if (data === false) {
    return <NoData className={clsx(className, styles.root)} classes={{placeholder: noDataStyle}}/>
  }

  if (showPrompt) {
    return <Alert
      severity="info"
      action={
        <Button
          color="inherit"
          size="small"
          onClick={e => { setShowPrompt(false); loadSystem(data, refViewer); setAccepted(true) }}
        >
            YES
        </Button>
      }
      data-testid={testID}
    >
      {`Visualization is by default disabled for systems with more than ${sizeLimit} atoms. Do you wish to enable visualization for this system with ${nAtoms} atoms?`}
    </Alert>
  }
  if (loading) {
    return <Placeholder
      variant="rect"
      className={clsx(styles.root, className)}
      data-testid={testID}
      classes={{placeholder: placeHolderStyle}}
    />
  }

  const menuItems = [
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
  ]
  if (data?.cell) {
    menuItems.push(...[
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
      </MenuItem>,
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
      </MenuItem>,
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
    ])
  }

  return <Floatable
    data-testid={testID}
    className={clsx(styles.root, className)}
    float={fullscreen}
    onFloat={toggleFullscreen}
  >
    <div className={styles.column}>
      {fullscreen && <Typography className={styles.title} variant="h6">Structure</Typography>}
      {loading
        ? <Placeholder
          variant="rect"
          className={styles.switchContainer}
          classes={{placeholder: styles.switchPlaceholder}}
        />
        : <div className={styles.canvas} ref={refCanvas}></div>
      }
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
          {m_path &&
            <Action tooltip='View data in the archive' onClick={() => { history.push(m_path) }}>
              <ViewList/>
            </Action>
          }
          <Action tooltip='Options' onClick={openMenu}>
            <MoreVert/>
          </Action>
        </Actions>
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
          {menuItems}
        </Menu>
      </div>
    </div>
  </Floatable>
})

Structure.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  viewer: PropTypes.object, // Optional shared viewer instance.
  data: PropTypes.oneOfType([
    PropTypes.bool,
    PropTypes.object
  ]),
  options: PropTypes.object, // Viewer options
  materialType: PropTypes.oneOf(['bulk', '2D', '1D', 'molecule / cluster', 'unavailable']),
  structureType: PropTypes.oneOf(['original', 'conventional', 'primitive']),
  m_path: PropTypes.string, // Path of the structure data in the metainfo
  captureName: PropTypes.string, // Name of the file that the user can download
  positionsOnly: PropTypes.bool, // Whether to update only positions. This is much faster than loading the entire structure.
  sizeLimit: PropTypes.number, // Maximum system size before a prompt is shown
  bondLimit: PropTypes.number, // The size at which bonds are turned off
  placeHolderStyle: PropTypes.string, // The CSS class to apply for the Placeholder component.
  noDataStyle: PropTypes.string, // The CSS class to apply for the NoData component.
  /**
   * A RxJS Subject for efficient, non-persistent, position changes that bypass
   * rendering of the component. Should send messages that contain the new
   * atomic positions as a list.
  */
  positionsSubject: PropTypes.any,
  'data-testid': PropTypes.string
}
Structure.defaultProps = {
  captureName: 'structure',
  sizeLimit: 500,
  bondLimit: 50
}

export default withWebGLErrorHandler(withErrorHandler(Structure, 'Could not load structure.'))
