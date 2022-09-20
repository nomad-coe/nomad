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
import Floatable from './Floatable'
import NoData from './NoData'
import Placeholder from '../visualization/Placeholder'
import { Actions, Action } from '../Actions'
import { withErrorHandler, withWebGLErrorHandler } from '../ErrorHandler'
import { useHistory } from 'react-router-dom'
import { isEmpty, flattenDeep, isNil, merge, cloneDeep } from 'lodash'
import { Quantity } from '../../units'
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
  captureName,
  sizeLimit,
  bondLimit,
  disableLegend,
  selection,
  'data-testid': testID,
  placeHolderStyle,
  noDataStyle}
) => {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showLatticeConstants, setShowLatticeConstants] = useState(true)
  const [showCell, setShowCell] = useState(true)
  const [center, setCenter] = useState('COP')
  const [fit, setFit] = useState('full')
  const [wrap, setWrap] = useState(true)
  const [showPrompt, setShowPrompt] = useState(false)
  const [accepted, setAccepted] = useState(false)
  const [showBonds, setShowBonds] = useState(true)
  const [species, setSpecies] = useState()
  const [loading, setLoading] = useState(true)
  const throwError = useAsyncError()

  // Variables
  const history = useHistory()
  const open = Boolean(anchorEl)
  const refViewer = useRef(null)
  const refCanvasDiv = useRef(null)
  const originalCenter = useRef()
  const firstLoad = useRef(true)
  const styles = useStyles(classes)
  const hasSelection = !isEmpty(selection?.selection)
  const nAtoms = useMemo(() => data?.positions.length, [data])
  const hasCell = useMemo(() => {
    return !data?.cell
      ? false
      : !flattenDeep(data.cell).every(v => v === 0)
  }, [data])

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
    // Dependency on 'fit' will make this function run every time it changes
    // (useCallback should not really run the function when the dependency
    // changes, but something weird is going on here.). An extra check is made
    // to check for node change before setting up the canvas to the new div.
    if (node !== refCanvasDiv.current) {
      refViewer.current.changeHostElement(node)
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
          enabled: showBonds,
          smoothness: isHighQuality ? 145 : 130,
          outline: {enabled: isHighQuality}
        }
      }
      refViewer.current = new StructureViewer(undefined, options)

      // If there is no valid cell, the property is set as undefined
      if (!hasCell) {
        system.cell = undefined
      }

      // Determine the orientation and view centering based on material type and
      // the structure type.
      let centerValue = 'COP'
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
      if (structureType === 'conventional' || structureType === 'primitive') {
        if (materialType === 'bulk') {
          centerValue = 'COC'
          alignments = [
            ['up', 'c'],
            ['right', 'b']
          ]
        } else if (materialType === '2D') {
          centerValue = 'COC'
          alignments = [
            ['right', 'a'],
            ['up', 'b']
          ]
          rotations = [
            [1, 0, 0, -60]
          ]
        } else if (materialType === '1D') {
          centerValue = 'COC'
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

      refViewer.current.load(system)
      refViewer.current.setRotation([1, 0, 0, 0])
      refViewer.current.atoms()
      refViewer.current.bonds({enabled: showBonds})
      if (hasCell) {
        refViewer.current.wrap(wrap)
        refViewer.current.cell({enabled: showCell})
        refViewer.current.latticeConstants({enabled: showLatticeConstants})
        refViewer.current.align(alignments)
      }
      refViewer.current.rotate(rotations)
      refViewer.current.controls({resetOnDoubleClick: false})
      originalCenter.current = centerValue
      const fitValue = 'full'
      refViewer.current.center(centerValue)
      refViewer.current.fit(fitValue, fitMargin)
      setCenter(centerValue)
      setFit(fitValue)
      refViewer.current.render()
      refViewer.current.saveCameraReset()

      // Get a list of all species ordered by atomic number
      const species = Object.entries(refViewer.current.elements)
        .map(([label, info]) => ({
          label: label,
          color: info[0],
          radius: info[1],
          atomicNumber: refViewer.current.elementNumbers[label]
        }))
        .sort((a, b) => a.atomicNumber - b.atomicNumber)
      setSpecies(species)
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [structureType, materialType, hasCell])

  // Called whenever the system changes. Loads the structure asynchronously.
  useEffect(() => {
    if (!data) return
    setLoading(true)

    if (!accepted) {
      if (nAtoms > bondLimit) {
        setShowBonds(false)
      }
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
        setLoading(false)
        firstLoad.current = false
      })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data, nAtoms, accepted, sizeLimit, bondLimit])

  // Handles selections. If a selection is given, the selected atoms will be
  // higlighted. Additionally the view will be centered on the selection if this
  // is requested.
  useEffect(() => {
    if (firstLoad.current) {
      return
    }

    // No selection: show all
    refViewer.current.resetCamera()
    let center
    let fit
    if (isEmpty(selection?.selection)) {
      center = originalCenter.current
      fit = 'full'
      refViewer.current.center(center)
      refViewer.current.fit(fit, fitMargin)
      refViewer.current.atoms()
      refViewer.current.bonds({enabled: showBonds})
      refViewer.current.cell({enabled: showCell})
      refViewer.current.latticeConstants({enabled: showLatticeConstants})
    // Focused selection: show only selected atoms
    } else {
      center = selection.isSubsystem ? selection.focus : originalCenter.current
      fit = selection.isSubsystem ? selection.focus : 'full'
      refViewer.current.center(center)
      refViewer.current.fit(fit, fitMargin)
      refViewer.current.atoms([
        {opacity: 0},
        {opacity: 0.1, include: selection.transparent},
        {opacity: 1, include: selection.selection}
      ])
      refViewer.current.bonds({enabled: false})
      refViewer.current.cell({enabled: !selection.isSubsystem})
      refViewer.current.latticeConstants({enabled: !selection.isSubsystem})
      setFit(fit)
      setCenter(center)
    }
    refViewer.current.render()
    refViewer.current.saveCameraReset()
    setFit(fit)
    setCenter(center)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selection])

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
          onClick={e => {
            setShowPrompt(false)
            setAccepted(true)
          }}
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
        : <div className={styles.canvas} ref={refCanvas}>
          {!disableLegend &&
            <Species species={species} className={styles.legend}/>
          }
        </div>
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
        <MenuItem key='show-bonds'>
          <FormControlLabel
            control={
              <Checkbox
                checked={hasSelection ? false : showBonds}
                disabled={hasSelection}
                onChange={(event) => {
                  setShowBonds(!showBonds)
                  refViewer.current.bonds({enabled: !showBonds})
                  refViewer.current.render()
                }}
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
                onChange={(event) => {
                  setShowLatticeConstants(!showLatticeConstants)
                  refViewer.current.latticeConstants({enabled: !showLatticeConstants})
                  refViewer.current.render()
                }}
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
                onChange={(event) => {
                  setShowCell(!showCell)
                  refViewer.current.cell({enabled: !showCell})
                  refViewer.current.render()
                }}
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
                onChange={(event) => {
                  setWrap(!wrap)
                  refViewer.current.wrap(!wrap)
                  refViewer.current.bonds({enabled: showBonds})
                  refViewer.current.render()
                }}
                color='primary'
              />
            }
            label='Wrap positions'
          />
        </MenuItem>
        </Menu>
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
  options: PropTypes.object, // Viewer options
  materialType: PropTypes.oneOf(['bulk', '2D', '1D', 'molecule / cluster', 'unavailable']),
  structureType: PropTypes.oneOf(['original', 'conventional', 'primitive']),
  m_path: PropTypes.string, // Path of the structure data in the metainfo
  captureName: PropTypes.string, // Name of the file that the user can download
  sizeLimit: PropTypes.number, // Maximum system size before a prompt is shown
  bondLimit: PropTypes.number, // The size at which bonds are turned off
  disableLegend: PropTypes.bool, // Disable the legend showing the species list
  placeHolderStyle: PropTypes.string, // The CSS class to apply for the Placeholder component.
  noDataStyle: PropTypes.string, // The CSS class to apply for the NoData component.
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

function getRadius(radius) {
  return `${Math.sqrt(radius) * 12}px`
}

/**
 * Shows a list of atomic species.
 */
const useSpeciesStyles = makeStyles((theme) => {
  return {
    root: {
      padding: theme.spacing(0.5),
      overflowY: 'auto',
      overflowX: 'hidden',
      direction: 'rtl'
    },
    grid: {
      display: 'grid',
      gridTemplateColumns: 'min-content auto',
      gridGap: `${theme.spacing(0.1)}px ${theme.spacing(0.5)}px`,
      direction: 'ltr'
    },
    label: {
      justifySelf: 'start',
      alignSelf: 'center',
      textShadow: '0px 0px 5px #ffffff'
    },
    circle: {
      borderRadius: '100%',
      border: '1px solid #555',
      justifySelf: 'center',
      alignSelf: 'center'
    }
  }
})
export const Species = React.memo(({species, className}) => {
  const styles = useSpeciesStyles()
  return <div className={clsx(styles.root, className)}>
    <div className={styles.grid}>
      {species && species.map(item => {
        const radius = getRadius(item.radius)
        return <React.Fragment key={item.label}>
          <div
            className={styles.circle}
            style={{
              backgroundColor: item.color,
              height: radius,
              width: radius
            }}
          />
          <Typography className={styles.label}>{item.label}</Typography>
        </React.Fragment>
      })}
    </div>
  </div>
})

Species.propTypes = {
  /**
   * A RxJS Subject for efficient, non-persistent, position changes that bypass
   * rendering of the component. Should send messages that contain the new
   * atomic positions as a list.
  */
  species: PropTypes.arrayOf(PropTypes.shape({
    label: PropTypes.string.isRequired, // The label to show
    radius: PropTypes.number.isRequired, // Radius in Ångstrom
    color: PropTypes.string.isRequired // CSS color
  })),
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

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

/**
 * Retrieves the topology for a calculation given the index and archive.
 *
 * @param {object} index
 * @param {object} archive
 *
 * @return {undefined|object} If the given structure cannot be converted,
 * returns an empty object.
 */
export function getTopology(index, archive) {
  // If topology is explicitly stored, use it.
  let topology
  if (index?.results?.material?.topology) {
    topology = merge(
      cloneDeep(index?.results?.material?.topology),
      archive?.results?.material?.topology
    )
  }

  // Create topology map
  let id = 0
  const topologyMap = Object.fromEntries(topology.map(top => {
    const node_id = `/results/material/topology/${id++}`
    return [node_id, top]
  }))

  // Create topology tree by finding the root node and then recursively
  // replacing its children with the actual child instances.
  const root = topology.find(top => isNil(top.parent_system))
  const traverse = (node) => {
    if (!isEmpty(node?.child_systems)) {
      node.child_systems = node.child_systems.map(id => topologyMap[id])
      node.child_systems.forEach(child => traverse(child))
    }
  }
  traverse(root)

  return [root, topologyMap]
}
