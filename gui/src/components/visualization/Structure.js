import React, { useState, useEffect, useRef, useCallback } from 'react'
import PropTypes from 'prop-types'

import { StructureViewer } from '@lauri-codes/materia'
import { makeStyles } from '@material-ui/core/styles'
import IconButton from '@material-ui/core/IconButton'
import Tooltip from '@material-ui/core/Tooltip'
import Menu from '@material-ui/core/Menu'
import { Button, Dialog, DialogTitle, DialogContent, DialogActions } from '@material-ui/core'
import MenuItem from '@material-ui/core/MenuItem'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import FullscreenIcon from '@material-ui/icons/Fullscreen'
import FullscreenExitIcon from '@material-ui/icons/FullscreenExit'
import CameraAltIcon from '@material-ui/icons/CameraAlt'
import ReplayIcon from '@material-ui/icons/Replay'
import FormControlLabel from '@material-ui/core/FormControlLabel'
import Checkbox from '@material-ui/core/Checkbox'
import Box from '@material-ui/core/Box'
import { convert } from '../../utils'

export default function Structure(props) {
  const [loaded, setLoaded] = React.useState(props.viewer === undefined ? false : props.viewer.structure !== undefined)
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showBonds, setShowBonds] = useState(true)
  const [showLatticeConstants, setShowLatticeConstants] = useState(true)
  const [showCell, setShowCell] = useState(true)
  const [error, setError] = useState(null)
  const open = Boolean(anchorEl)

  const viewer = useRef(null)
  const refCanvas = useRef(null)

  const useStyles = makeStyles((theme) => {
    return {
      root: {
        width: '100%',
        height: 0,
        paddingBottom: `${(1 / props.aspectRatio) * 100}%`, // CSS hack for fixed aspect ratio
        position: 'relative'
      },
      rootInner: {
        position: 'absolute',
        top: 0,
        right: 0,
        bottom: 0,
        left: 0
      },
      container: {
        display: 'flex',
        width: '100%',
        height: '100%',
        flexDirection: 'column',
        backgroundColor: 'white'
      },
      header: {
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
        marginBottom: theme.spacing(2),
        display: error === null ? 'block' : 'none'
      },
      errorContainer: {
        flex: 1,
        zIndex: 0,
        minHeight: 0, // added min-height: 0 to allow the item to shrink to fit inside the container.
        marginBottom: theme.spacing(2),
        alignItems: 'center',
        justifyContent: 'center',
        display: error === null ? 'none' : 'flex'
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
  const classes = useStyles(props)

  // In order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const measuredRef = useCallback(node => {
    refCanvas.current = node
    if (node === null) {
      return
    }
    if (viewer.current === null) {
      return
    }
    viewer.current.changeHostElement(node, true, true)
  }, [])

  // Run only on first render to initialize the viewer. See the viewer
  // documentation for details on the meaning of different options:
  // https://lauri-codes.github.io/materia/viewers/structureviewer
  useEffect(() => {
    let options
    if (props.options === undefined) {
      options = {
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
      options = props.options
    }
    if (props.viewer === undefined) {
      viewer.current = new StructureViewer(undefined, options)
    } else {
      viewer.current = props.viewer
      viewer.current.setOptions(options, false, false)
    }
    if (refCanvas.current !== null) {
      viewer.current.changeHostElement(refCanvas.current, false, false)
    }
  }, [])

  // For loading a structure
  useEffect(() => {
    // Check data validity. Display error message if data is invalid or too big
    // to display.
    if (props.system === undefined) {
      return
    }
    let nAtoms = props.system.atom_species.length
    if (nAtoms >= props.sizeLimit) {
      setError('Visualization is disabled due to large system size.')
      return
    }

    // If a structure has already been loaded and a viewer has been pre-created,
    // it is assumed that only the positions need to be updated on data change.
    if (props.viewer !== undefined && loaded) {
      let positions = convert(props.system.atom_positions, 'm', 'angstrom')
      viewer.current.setPositions(positions)
      return
    }

    // Systems with cell are centered on the cell center and oriented is defined
    // by the cell vectors.
    let cell = props.system.lattice_vectors
    if (cell !== undefined) {
      viewer.current.setOptions({layout: {
        viewCenter: 'COC',
        viewRotation: {
          align: {
            top: 'c',
            right: 'b'
          },
          rotations: [
            [0, 1, 0, 60],
            [1, 0, 0, 30]
          ]
        }
      }})
    // Systems without cell are centered on the center of positions
    } else {
      viewer.current.setOptions({layout: {
        viewCenter: 'COP',
        viewRotation: {
          rotations: [
            [0, 1, 0, 60],
            [1, 0, 0, 30]
          ]
        }
      }})
    }
    var system = {
      'atomicNumbers': props.system.atom_species,
      'cell': convert(props.system.lattice_vectors, 'm', 'angstrom'),
      'positions': convert(props.system.atom_positions, 'm', 'angstrom'),
      'pbc': props.system.configuration_periodic_dimensions
    }
    viewer.current.load(system)
    viewer.current.fitToCanvas()
    viewer.current.saveReset()
    viewer.current.reset()
    setLoaded(true)
  }, [])

  // Viewer settings
  useEffect(() => {
    viewer.current.setOptions({bonds: {enabled: showBonds}})
  }, [showBonds])

  useEffect(() => {
    viewer.current.setOptions({latticeConstants: {enabled: showLatticeConstants}})
  }, [showLatticeConstants])

  useEffect(() => {
    viewer.current.setOptions({cell: {enabled: showCell}})
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
    viewer.current.takeScreenShot(props.captureName)
  }, [props.captureName])

  const handleReset = useCallback(() => {
    viewer.current.reset()
    viewer.current.fitToCanvas()
    viewer.current.render()
  }, [])

  let content = <Box className={classes.container}>
    <div className={classes.header}>
      <div className={classes.spacer}></div>
      <Tooltip title="Reset view">
        <IconButton className={classes.iconButton} onClick={handleReset} disabled={error}>
          <ReplayIcon />
        </IconButton>
      </Tooltip>
      <Tooltip
        title="Toggle fullscreen">
        <IconButton className={classes.iconButton} onClick={toggleFullscreen} disabled={error}>
          {fullscreen ? <FullscreenExitIcon /> : <FullscreenIcon />}
        </IconButton>
      </Tooltip>
      <Tooltip title="Capture image">
        <IconButton className={classes.iconButton} onClick={takeScreencapture} disabled={error}>
          <CameraAltIcon />
        </IconButton>
      </Tooltip>
      <Tooltip title="Options">
        <IconButton className={classes.iconButton} onClick={openMenu} disabled={error}>
          <MoreVertIcon />
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
    <div className={classes.viewerCanvas} ref={measuredRef}></div>
    <div className={classes.errorContainer}><div className={classes.errorMessage}>{error}</div></div>
    }
  </Box>

  return (
    <Box className={classes.root}>
      <Box className={classes.rootInner}>
        {fullscreen ? '' : content}
      </Box>
      <Dialog maxWidth="lg" fullWidth open={fullscreen}>
        <DialogTitle>{'Structure'}</DialogTitle>
        <DialogContent>
          <Box className={classes.root}>
            <Box className={classes.rootInner}>
              {fullscreen ? content : ''}
            </Box>
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setFullscreen(false)}>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </Box>
  )
}

Structure.propTypes = {
  system: PropTypes.object,
  options: PropTypes.object,
  captureName: PropTypes.string,
  viewer: PropTypes.any,
  aspectRatio: PropTypes.number,
  sizeLimit: PropTypes.number
}
Structure.defaultProps = {
  aspectRatio: 4 / 3,
  captureName: 'structure',
  sizeLimit: 300
}
