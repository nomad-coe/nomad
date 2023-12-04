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
t* See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React, { useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import { useHistory } from 'react-router-dom'
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
  ViewList,
  GetApp
} from '@material-ui/icons'
import { DownloadSystemMenu, WrapModeRadio } from '../buttons/DownloadSystemButton'
import Floatable from './Floatable'
import NoData from './NoData'
import Placeholder from './Placeholder'
import { Actions, Action } from '../Actions'
import { withErrorHandler, withWebGLErrorHandler } from '../ErrorHandler'
import { isNil } from 'lodash'
import { Quantity } from '../units/Quantity'

/**
 * Used to control a 3D system visualization that is implemented in the
 * 'children' prop. This allows for an easier change of visualization
 * implementation without changing the control layout.
 */
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
    title: {
      marginBottom: theme.spacing(1)
    },
    menuItem: {
      margin: theme.spacing(2),
      marginTop: theme.spacing(1),
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
const StructureBase = React.memo(({
  wrapMode,
  onWrapModeChange,
  disableWrapMode,
  showLatticeConstants,
  onShowLatticeConstants,
  disableShowLatticeConstants,
  showBonds,
  onShowBonds,
  disableShowBonds,
  showCell,
  onShowCell,
  disableShowCell,
  onAccept,
  onReset,
  onTakeScreenshot,
  onFullscreen,
  float,
  species,
  loading,
  prompt,
  noData,
  className,
  classes,
  captureName,
  disableLegend,
  disableFileDownload,
  entryId,
  path,
  m_path,
  'data-testid': testID,
  children}
) => {
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [systemAnchorEl, setSystemAnchorEl] = React.useState(null)
  const history = useHistory()
  const open = Boolean(anchorEl)
  const styles = useStyles(classes)
  const downloadDisabled = isNil(entryId) || isNil(path) || disableFileDownload

  const setShowBonds = useCallback((value, render = false) => {
    onShowBonds && onShowBonds(value, render)
  }, [onShowBonds])

  const handleWrapModeChange = useCallback((value, showBonds, render = false) => {
    onWrapModeChange && onWrapModeChange(value, showBonds, render)
  }, [onWrapModeChange])

  const handleShowLatticeConstants = useCallback((value, render = false) => {
    onShowLatticeConstants && onShowLatticeConstants(value, render)
  }, [onShowLatticeConstants])

  const setShowCell = useCallback((value, render = false) => {
    onShowCell && onShowCell(value, render)
  }, [onShowCell])

  const handleFullscreen = useCallback(() => {
    onFullscreen && onFullscreen()
  }, [onFullscreen])

  const handleReset = useCallback(() => {
    onReset && onReset()
  }, [onReset])

  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])

  const openSystemMenu = useCallback((event) => {
    setSystemAnchorEl(event.currentTarget)
  }, [])

  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

  const takeScreencapture = useCallback(() => {
    onTakeScreenshot && onTakeScreenshot(captureName)
  }, [onTakeScreenshot, captureName])

  // Decides the possible overlay to display
  const overlay = useMemo(() => {
    if (prompt) {
      return <Alert
        severity="info"
        action={
          <Button color="inherit" size="small" onClick={() => onAccept()}>
              YES
          </Button>
        }
      >
        {prompt}
      </Alert>
    } else if (loading) {
      return <Placeholder margin={0}/>
    } else if (noData) {
      return <NoData margin={0}/>
    }
  }, [loading, prompt, noData, onAccept])

  return <Floatable
    data-testid={testID}
    className={clsx(styles.root, className)}
    float={float}
  >
    <div className={styles.container}>
      <div className={styles.content}>
        {overlay}
      </div>
      {/* This contains the actual visualization. Needs to be always in the DOM so we can
      load the structure into it properly. */}
      <div
        className={styles.content}
        style={{visibility: overlay ? 'hidden' : 'visible'}}
      >
        {float && <Typography className={styles.title} variant="h6">Structure</Typography>}
        <div className={styles.canvas}>
          {children}
          {!disableLegend &&
            <Species species={species} className={styles.legend}/>
          }
        </div>
        <div className={styles.header}>
          <Actions>
            <Action tooltip='Reset view' onClick={handleReset}>
              <Replay/>
            </Action>
            <Action tooltip='Toggle fullscreen' onClick={handleFullscreen}>
              {float ? <FullscreenExit/> : <Fullscreen/>}
            </Action>
            <Action
              tooltip={downloadDisabled
                ? 'File download not available for this system'
                : 'Download system geometry as a file'
              }
              onClick={openSystemMenu}
              disabled={downloadDisabled}
            >
              <GetApp/>
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
          <DownloadSystemMenu
            entryId={entryId}
            path={path}
            anchorEl={systemAnchorEl}
            onClose={() => setSystemAnchorEl(null)}
          />
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
                    disabled={disableShowBonds}
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
                    checked={showLatticeConstants}
                    disabled={disableShowLatticeConstants}
                    onChange={(event) => { handleShowLatticeConstants(!showLatticeConstants, true) }}
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
                    disabled={disableShowCell}
                    onChange={(event) => { setShowCell(!showCell, true) }}
                    color='primary'
                  />
                }
                label='Show simulation cell'
              />
            </MenuItem>
            <WrapModeRadio value={wrapMode} onChange={handleWrapModeChange} disabled={disableWrapMode} className={styles.menuItem}/>
          </Menu>
        </div>
      </div>
    </div>
  </Floatable>
})

StructureBase.propTypes = {
  wrapMode: PropTypes.bool,
  onWrapModeChange: PropTypes.func,
  disableWrapMode: PropTypes.bool,
  showLatticeConstants: PropTypes.bool,
  onShowLatticeConstants: PropTypes.func,
  disableShowLatticeConstants: PropTypes.bool,
  showBonds: PropTypes.bool,
  onShowBonds: PropTypes.func,
  disableShowBonds: PropTypes.bool,
  showCell: PropTypes.bool,
  onShowCell: PropTypes.func,
  disableShowCell: PropTypes.bool,
  prompt: PropTypes.string, // Prompt shown before showing the visualization
  onAccept: PropTypes.func, // Callback for answering yes to the prompt
  onReset: PropTypes.func,
  onTakeScreenshot: PropTypes.func,
  onFullscreen: PropTypes.func,
  species: PropTypes.array,
  loading: PropTypes.bool,
  float: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  noData: PropTypes.bool,
  canvasID: PropTypes.string,
  disableFileDownload: PropTypes.bool, // Used to disable the file download button
  entryId: PropTypes.string,
  path: PropTypes.string,
  m_path: PropTypes.string, // Path of the structure data in the metainfo
  captureName: PropTypes.string, // Name of the file that the user can download
  disableLegend: PropTypes.bool, // Disable the legend showing the species list
  'data-testid': PropTypes.string,
  children: PropTypes.node
}
StructureBase.defaultProps = {
  captureName: 'system'
}

export default withWebGLErrorHandler(withErrorHandler('Could not load structure.')(StructureBase))

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
  const normalized = useMemo(() => {
    let normalized
    try {
      normalized = species
        .map(x => ({...x, radius: (x.radius)}))
        .sort((a, b) => (a.atomicNumber - b.atomicNumber))
    } catch {}
    return normalized
  }, [species])

  return normalized
    ? <div className={clsx(styles.root, className)}>
      <div className={styles.grid}>
        {species && normalized.map(item => {
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
    : null
})

Species.propTypes = {
  species: PropTypes.arrayOf(PropTypes.shape({
    label: PropTypes.string.isRequired, // The label to show
    radius: PropTypes.number.isRequired, // Radius in Ångstrom
    color: PropTypes.string.isRequired, // CSS color
    atomicNumber: PropTypes.number.isRequired // Atomic number for sorting
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
