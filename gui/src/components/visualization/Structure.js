import React, { useState, useEffect, useRef, useCallback } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import {
  Box,
  Checkbox,
  Menu,
  MenuItem,
  Button,
  IconButton,
  Tooltip,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
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

export default function Structure(props) {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showBonds, setShowBonds] = useState(true)
  const [showLatticeConstants, setShowLatticeConstants] = useState(true)
  const [showCell, setShowCell] = useState(true)
  const [error, setError] = useState(null)

  // Variables
  const open = Boolean(anchorEl)
  const viewer = useRef(null)
  const refCanvas = useRef(null)

  // Styles
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
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Called only on first render to load the given structure.
  useEffect(() => {
    if (props.system === undefined) {
      return
    }

    if (props.positionsOnly) {
      viewer.current.setPositions(props.system.positions)
      return
    }

    let nAtoms = props.system.species.length
    if (nAtoms >= props.sizeLimit) {
      setError('Visualization is disabled due to large system size.')
      return
    }

    // Systems with cell are centered on the cell center and orientation is defined
    // by the cell vectors.
    let cell = props.system.cell
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
    viewer.current.load(props.system)
    viewer.current.fitToCanvas()
    viewer.current.saveReset()
    viewer.current.reset()
  // eslint-disable-next-line react-hooks/exhaustive-deps
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
          <Replay />
        </IconButton>
      </Tooltip>
      <Tooltip
        title="Toggle fullscreen">
        <IconButton className={classes.iconButton} onClick={toggleFullscreen} disabled={error}>
          {fullscreen ? <FullscreenExit /> : <Fullscreen />}
        </IconButton>
      </Tooltip>
      <Tooltip title="Capture image">
        <IconButton className={classes.iconButton} onClick={takeScreencapture} disabled={error}>
          <CameraAlt />
        </IconButton>
      </Tooltip>
      <Tooltip title="Options">
        <IconButton className={classes.iconButton} onClick={openMenu} disabled={error}>
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
    <div className={classes.viewerCanvas} ref={measuredRef}></div>
    <div className={classes.errorContainer}><div className={classes.errorMessage}>{error}</div></div>
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
  viewer: PropTypes.object, // Optional shared viewer instance.
  id: PropTypes.string, // Id for the visualized structure.
  system: PropTypes.object, // The system to display as section_system
  options: PropTypes.object, // Viewer options
  captureName: PropTypes.string, // Name of the file that the user can download
  aspectRatio: PropTypes.number, // Fixed aspect ratio for the viewer canvas
  sizeLimit: PropTypes.number, // Maximum number of atoms to attempt to display
  positionsOnly: PropTypes.bool // Whether to update only positions. This is much faster than loading the entire structure.
}
Structure.defaultProps = {
  aspectRatio: 4 / 3,
  captureName: 'structure',
  sizeLimit: 300
}
