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
import React, { useState, useMemo, useEffect, useRef, useCallback } from 'react'
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
import { Species } from './StructureBase'
import Floatable from './Floatable'
import NoData from './NoData'
import Placeholder from '../visualization/Placeholder'
import { Actions, Action } from '../Actions'
import { withErrorHandler, withWebGLErrorHandler } from '../ErrorHandler'
import { useHistory } from 'react-router-dom'
import { isEmpty, flattenDeep } from 'lodash'
import { Quantity } from '../units/Quantity'
import { delay } from '../../utils'
import { useAsyncError } from '../../hooks'
import clsx from 'clsx'

/**
 * Used to show atomistic systems in an interactive 3D viewer based on the
 * 'materia'-library.
 */
const fitMargin = 0.5
const useStyles = makeStyles((theme) => {
  return {
    root: {},
    container: {
      position: 'relative',
      width: '100%',
      height: '100%'
    },
    content: {
      position: 'absolute',
      left: 0,
      right: 0,
      top: 0,
      bottom: 0,
      display: 'flex',
      flexDirection: 'column'
    },
    header: {
      display: 'flex',
      flexDirection: 'row',
      zIndex: 1
    },
    legend: {
      zIndex: 10,
      position: 'absolute',
      left: 0,
      top: 0,
      bottom: 0
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
      position: 'relative',
      flex: 1,
      zIndex: 0,
      minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
      marginBottom: theme.spacing(1)
    }
  }
})
const Structure = React.memo(({
  className,
  classes,
  data,
  structuralType,
  cellType,
  m_path,
  captureName,
  sizeLimit,
  bondLimit,
  disableLegend,
  selection,
  'data-testid': testID
}) => {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showLatticeConstants, setShowLatticeConstantsState] = useState(true)
  const [showCell, setShowCellState] = useState(true)
  const [center, setCenter] = useState('COP')
  const [fit, setFit] = useState('full')
  const [wrap, setWrapState] = useState(true)
  const [showPrompt, setShowPrompt] = useState(false)
  const [accepted, setAccepted] = useState(false)
  const [showBonds, setShowBondsState] = useState(true)
  const [species, setSpecies] = useState()
  const [loading, setLoading] = useState(true)
  const [ready, setReady] = useState(false)
  const throwError = useAsyncError()
  const noData = data === false

  const setShowBonds = useCallback((value, render = false) => {
    setShowBondsState(value)
    if (refViewer?.current) {
      try {
        refViewer.current.bonds({enabled: value})
        if (render) refViewer.current.render()
      } catch (e) {
      }
    }
  }, [])

  const setWrap = useCallback((value, showBonds, render = false) => {
    setWrapState(value)
    if (refViewer?.current) {
      try {
          refViewer.current.wrap(value)
          refViewer.current.bonds({enabled: showBonds})
        if (render) refViewer.current.render()
      } catch (e) {
      }
    }
  }, [])

  const setShowLatticeConstants = useCallback((value, render = false) => {
    setShowLatticeConstantsState(value)
    if (refViewer?.current) {
      try {
          refViewer.current.latticeConstants({enabled: value})
        if (render) refViewer.current.render()
      } catch (e) {
      }
    }
  }, [])

  const setShowCell = useCallback((value, render = false) => {
    setShowCellState(value)
    if (refViewer?.current) {
      try {
          refViewer.current.cell({enabled: value})
        if (render) refViewer.current.render()
      } catch (e) {
      }
    }
  }, [])

  // Variables
  const history = useHistory()
  const open = Boolean(anchorEl)
  const readyRef = useRef(true)
  const refViewer = useRef(null)
  const refCanvasDiv = useRef(null)
  const styles = useStyles(classes)
  const hasSelection = !isEmpty(selection?.selection)
  const nAtoms = useMemo(() => data?.positions?.length, [data])
  const hasCell = useMemo(() => {
    const hasCell = !data?.cell
      ? false
      : !flattenDeep(data.cell).every(v => v === 0)
    // If there is no valid cell, the property is set as undefined
    if (!hasCell && data) {
      data.cell = undefined
    }
    return hasCell
  }, [data])

  // In order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const refCanvas = useCallback(node => {
    // Dependency on 'fit' will make this function run every time it changes
    // (useCallback should not really run the function when the dependency
    // changes, but something weird is going on here.). An extra check is made
    // to check for node change before setting up the canvas to the new div.
    if (node !== null && refViewer.current !== null && node !== refCanvasDiv.current) {
      refViewer.current.changeHostElement(node, true)
      refViewer.current.fit(fit, fitMargin)
      refViewer.current.render()
    }
    refCanvasDiv.current = node
  }, [fit])

  /**
   * Used to asynchronously load a new viewer and system. Since loading the
   * structure is a long-taking operation, running this function later in the
   * event queue allows the component to perform state updates (e.g. loading
   * placeholder) while the viewer is loading.
   */
  const loadSystem = useCallback(async (system, refViewer) => {
    await delay(() => {
      // Initialize the viewer. A new viewer is used for each system as the
      // render settings may differ.
      const nAtoms = system?.positions?.length
      const isHighQuality = nAtoms <= 3000
      const options = {
        renderer: {
          backgroundColor: ['#ffffff', 1],
          antialias: {enabled: isHighQuality},
          pixelRatioScale: 1,
          shadows: {enabled: false}
        },
        atoms: {
          smoothness: isHighQuality ? 165 : 150,
          outline: {enabled: isHighQuality}
        },
        bonds: {
          enabled: false,
          smoothness: isHighQuality ? 145 : 130,
          outline: {enabled: isHighQuality}
        }
      }
      refViewer.current = new StructureViewer(refCanvasDiv.current, options)
      refViewer.current.load(system)
      refViewer.current.controls({resetOnDoubleClick: false})

      // Get a list of all species ordered by atomic number
      const species = Object.entries(refViewer.current.elements)
        .map(([label, info]) => ({
          label: label,
          color: info[0],
          radius: info[1],
          atomicNumber: refViewer.current.elementNumbers[label]
        }))
      setSpecies(species)
    })
  }, [])

  // Called whenever the system changes. Loads the structure asynchronously.
  useEffect(() => {
    if (!data) return
    readyRef.current = false
    setReady(false)
    setLoading(true)

    // If the system is very large, ask the user first for permission to attempt
    // to visualize it.
    if (!accepted) {
      if (nAtoms > sizeLimit) {
        setShowPrompt(true)
        return
      }
    }

    // Remember to catch since react error boundaries do not automatically catch
    // from async calls.
    loadSystem(data, refViewer)
      .catch(throwError)
      .finally(() => {
        setReady(true)
        readyRef.current = true
      })
  }, [data, nAtoms, accepted, sizeLimit, throwError, loadSystem])

  // Once the system is loaded, this effect will determine the final visual
  // layout. By monitoring the ready-state, the selections are updated correctly
  // once the system is loaded. By also monitoring the readyRef-reference, this
  // effect can know the status within the same render cycle as well.
  useEffect(() => {
    if (!ready || !readyRef.current) {
      return
    }

    // Reset camera and system rotations
    refViewer.current.resetCamera()
    refViewer.current.setRotation([1, 0, 0, 0])

    // Determine the orientation and view centering based on material type and
    // the structure type.
    let center = 'COP'
    let fit = 'full'
    let showBondsValue = nAtoms < bondLimit
    let rotations = [
      [0, 1, 0, 60],
      [1, 0, 0, 30]
    ]
    let alignments = hasCell
      ? [
        ['up', 'c'],
        ['right', 'b']
      ]
      : undefined
    if (cellType === 'conventional' || cellType === 'primitive') {
      if (structuralType === 'bulk') {
        center = 'COC'
        alignments = [
          ['up', 'c'],
          ['right', 'b']
        ]
      } else if (structuralType === '2D') {
        center = 'COC'
        alignments = [
          ['right', 'a'],
          ['up', 'b']
        ]
        rotations = [
          [1, 0, 0, -60]
        ]
      } else if (structuralType === '1D') {
        center = 'COC'
        alignments = [
          ['right', 'a']
        ]
        rotations = [
          [1, 0, 0, -60],
          [0, 1, 0, 30],
          [1, 0, 0, 30]
        ]
      }
    }

    // No selection: show all
    if (isEmpty(selection?.selection)) {
      refViewer.current.atoms()
      showBondsValue = nAtoms < bondLimit
      setShowBonds(showBondsValue, false)
      if (hasCell) {
        setShowCell(true)
        setShowLatticeConstants(true)
      }
    // Focused selection: show only selected atoms
    } else {
      center = selection.isSubsystem ? selection.focus : center
      fit = selection.isSubsystem ? selection.focus : 'full'
      refViewer.current.atoms([
        {opacity: 0},
        {opacity: 0.1, include: selection.transparent},
        {opacity: 1, include: selection.selection}
      ])
      // Bonds are not shown for selections
      if (!selection.isSubsystem) {
        showBondsValue = false
        setShowBonds(showBondsValue, false)
      }
      if (hasCell) {
        setShowCell(!selection.isSubsystem)
        setShowLatticeConstants(!selection.isSubsystem)
      }
    }
    if (hasCell) {
      setWrap(true, showBondsValue, false)
      refViewer.current.align(alignments)
    }
    refViewer.current.rotate(rotations)
    refViewer.current.center(center)
    refViewer.current.fit(fit, fitMargin)
    refViewer.current.render()
    refViewer.current.saveCameraReset()
    setFit(fit)
    setCenter(center)
    setLoading(false)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selection, ready])

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
    refViewer.current.resetCamera()
    refViewer.current.center(center)
    refViewer.current.fit(fit, fitMargin)
    refViewer.current.render()
  }, [center, fit])

  // Decides the possible overlay to display
  const overlay = useMemo(() => {
    if (showPrompt) {
      return <Alert
        severity="info"
        action={
          <Button color="inherit" size="small" onClick={() => {
            setShowPrompt(false)
            setAccepted(true)
          }}>
              YES
          </Button>
        }
      >
        {`Visualization is by default disabled for systems with more than ${sizeLimit} atoms. Do you wish to enable visualization for this system with ${nAtoms} atoms?`}
      </Alert>
    } else if (loading) {
      return <Placeholder margin={0}/>
    } else if (noData) {
      return <NoData margin={0}/>
    }
    return null
  }, [loading, nAtoms, noData, showPrompt, sizeLimit])

  return <Floatable
    data-testid={testID}
    className={clsx(styles.root, className)}
    float={fullscreen}
    onFloat={toggleFullscreen}
  >
    <div className={styles.container}>
      {/* Overlay that is displayed on top of the viewer */}
      <div className={styles.content}>
        {overlay}
      </div>
      {/* This contains the actual visualization. Needs to be always in the DOM so we can
      load the structure into it properly. */}
      <div
        className={styles.content}
        style={{visibility: overlay ? 'hidden' : 'visible'}}
      >
        {fullscreen && <Typography className={styles.title} variant="h6">Structure</Typography>}
        <div className={styles.canvas} ref={refCanvas}>
          {!disableLegend &&
            <Species species={species} className={styles.legend}/>
          }
        </div>
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
          <MenuItem key='show-bonds'>
            <FormControlLabel
              control={
                <Checkbox
                  checked={hasSelection ? false : showBonds}
                  disabled={hasSelection}
                  onChange={(event) => { setShowBonds(!showBonds, true) }}
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
                  checked={!hasCell ? false : (hasSelection ? !selection.isSubsystem : showLatticeConstants)}
                  disabled={!hasCell ? true : hasSelection}
                  onChange={(event) => { setShowLatticeConstants(!showLatticeConstants, true) }}
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
                  checked={!hasCell ? false : (hasSelection ? !selection.isSubsystem : showCell)}
                  disabled={!hasCell ? true : hasSelection}
                  onChange={(event) => { setShowCell(!showCell, true) }}
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
                  checked={!hasCell ? false : wrap}
                  disabled={!hasCell ? true : hasSelection}
                  onChange={(event) => { setWrap(!wrap, showBonds, true) }}
                  color='primary'
                />
              }
              label='Wrap positions'
            />
          </MenuItem>
          </Menu>
        </div>
      </div>
    </div>
  </Floatable>
})

Structure.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  data: PropTypes.oneOfType([
    PropTypes.bool,
    PropTypes.object
  ]),
  structuralType: PropTypes.oneOf(['bulk', 'surface', '2D', '1D', 'molecule / cluster', 'unavailable']),
  cellType: PropTypes.oneOf(['original', 'conventional', 'primitive']),
  m_path: PropTypes.string, // Path of the structure data in the metainfo
  captureName: PropTypes.string, // Name of the file that the user can download
  sizeLimit: PropTypes.number, // Maximum system size before a prompt is shown
  bondLimit: PropTypes.number, // The size at which bonds are turned off
  disableLegend: PropTypes.bool, // Disable the legend showing the species list
  selection: PropTypes.shape({
    /**
     * The atom indices to show.
     */
    selection: PropTypes.arrayOf(PropTypes.number),
    /**
     * The atom indices to which the view is focused on.
     */
    focus: PropTypes.arrayOf(PropTypes.number),
    /**
     * The atom indices to display transparently.
     */
    transparent: PropTypes.arrayOf(PropTypes.number),
    /**
     * Is the selection visualizing a subsystem, e.g. a single molecule?
     */
    isSubsystem: PropTypes.bool
  }),
  'data-testid': PropTypes.string
}
Structure.defaultProps = {
  captureName: 'structure',
  sizeLimit: 500,
  bondLimit: 50
}

export default withWebGLErrorHandler(withErrorHandler('Could not load structure.')(Structure))

/**
 * Converts the given structure in the format used by 'results' into the format
 * used by the materia-library.
 *
 * @param {object} structure.
 *
 * @return {undefined|object} If the given structure cannot be converted,
 * returns an empty object.
 */
export function toMateriaStructure(structure) {
  if (!structure) {
    return undefined
  }

  try {
    // Resolve atomic species using the labels and their mapping to chemical
    // elements.
    const speciesMap = new Map(structure.species.map(s => [s.name, s.chemical_symbols[0]]))

    const structMateria = {
      species: structure.species_at_sites.map(x => speciesMap.get(x)),
      cell: structure.lattice_vectors
        ? new Quantity(structure.lattice_vectors, 'meter').to('angstrom').value()
        : undefined,
      positions: new Quantity(structure.cartesian_site_positions, 'meter').to('angstrom').value(),
      fractional: false,
      pbc: structure.dimension_types ? structure.dimension_types.map((x) => !!x) : undefined
    }
    return structMateria
  } catch (error) {
    return {}
  }
}
