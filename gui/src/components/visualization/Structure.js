import React, { useState, useEffect, useRef, useLayoutEffect } from 'react'
import PropTypes from 'prop-types'

import { StructureViewer } from '@lauri-codes/materia'
import { makeStyles } from '@material-ui/core/styles'
import { Paper } from '@material-ui/core'
import IconButton from '@material-ui/core/IconButton'
import Menu from '@material-ui/core/Menu'
import MenuItem from '@material-ui/core/MenuItem'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import FullscreenIcon from '@material-ui/icons/Fullscreen'
import CameraAltIcon from '@material-ui/icons/CameraAlt'
import ReplayIcon from '@material-ui/icons/Replay'
import FormControlLabel from '@material-ui/core/FormControlLabel'
import ClickAwayListener from '@material-ui/core/ClickAwayListener'
import Checkbox from '@material-ui/core/Checkbox'
import Box from '@material-ui/core/Box'

export default function Structure(props) {
  // States
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [fullscreen, setFullscreen] = useState(false)
  const [showBonds, setShowBonds] = useState(true)
  const [showAxis, setShowAxis] = useState(true)
  const [dimensions, setDimensions] = useState({ width: 1, height: 1 })
  const open = Boolean(anchorEl)

  // References
  const viewer = useRef(null)
  const refViewerCanvas = useRef(null)
  const refContainer = useRef(null)

  // Styles
  const useStyles = makeStyles((theme) => {
    // Calculate sizing
    let maxSize = 0.9
    let width, height
    let aspectRatio = dimensions.height / dimensions.width
    let windowRatio = window.innerHeight / window.innerWidth
    if (aspectRatio >= windowRatio) {
      height = maxSize * window.innerHeight
      width = height / aspectRatio
    } else {
      width = maxSize * window.innerWidth
      height = width * aspectRatio
    }

    return {
      root: {
        width: '100%',
        boxSizing: 'border-box' // This should be a global default?
      },
      container: {
        boxSizing: 'border-box', // This should be a global default?
        display: 'flex',
        width: '100%',
        height: '100%',
        flexDirection: 'column',
        backgroundColor: 'white'
      },
      fullscreen: {
        boxSizing: 'border-box', // This should be a global default?
        padding: theme.spacing(2),
        position: 'fixed',
        zIndex: 2000,
        left: '50%',
        top: '50%',
        width: width + 'px',
        height: height + 'px',
        transform: 'translate(-50%, -50%)'
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
        marginBottom: theme.spacing(2)
      },
      iconButton: {
        backgroundColor: 'white'
      }
    }
  })
  const classes = useStyles(props)

  // Event handlers
  const handleOpenMenu = (event) => {
    setAnchorEl(event.currentTarget)
  }

  const handleCloseMenu = () => {
    setAnchorEl(null)
  }

  const handleFullscreen = () => {
    setFullscreen(!fullscreen)
  }

  const handleScreencapture = () => {
    viewer.current.takeScreenShot(props.calcId)
  }

  const handleClickAway = (event) => {
    // For some reason clicking the screen capture button will cause the browser
    // to fire a click event during the file download. This case is detected by
    // checking the event target properties.
    let notDownload = event.target.download === undefined
    if (fullscreen && notDownload) {
      setFullscreen(false)
    }
  }

  const handleReset = () => {
    viewer.current.reset()
    viewer.current.fitToCanvas()
    viewer.current.render()
  }

  // Used to fetch the initial dimensions of the structure viewer. This fixes
  // the aspect ratio that will be used in fullscreen mode.
  useLayoutEffect(() => {
    if (refContainer.current) {
      setDimensions({
        width: refContainer.current.offsetWidth,
        height: refContainer.current.offsetHeight
      })
    }
  }, [])

  // Run only on first render to initialize the viewer
  useEffect(() => {
    let options = {
      view: {
        autoResize: true, // Automatically fit canvas to the root visualization element
        autoFit: true, // Automatically set the zoom to fit the structure into the canvas
        fitMargin: 0.0 // Margin between visualization and canvas boundary
      },
      controls: {
        enableZoom: true, // Enable zoom with mouse wheel/pinch
        enablePan: true, // Enable pan with mouse/touch
        enableRotate: true // Enable rotation
      },
      structure: {
        showParam: true, // Show lattice parameters
        showCopies: true, // Show periodic copies at cell boundary
        showShadows: false, // Enable shadows: requires a bit more from GPU
        allowRepeat: false, // Allow repeating the structure
        showCell: true, // Show unit cell wireframe
        wrap: false, // Wrap atoms to unit cell
        createLegend: false, // Show atom labels
        showLegend: true, // Show atom labels
        showOptions: false, // Show the options toolbar
        showBonds: true, // Show or hide bonds
        radiusScale: 0.8, // Scaling factor for atomic radii
        bondScale: 1.2, // Scaling factor for the automatic bond detection
        viewCenter: 'COC' // The rotation and view center position. Valid values are: "COP" (center of position), "COC" center of cell, or a custom position given as array of three cartesian coordinates.
      },
      style: {
        backgroundColor: [0xffffff, 1]
      }
    }
    viewer.current = new StructureViewer(refViewerCanvas.current, options)
    var bulk = {
      'atomicNumbers': [11, 17, 11, 17, 11, 17, 11, 17],
      'cell': [
        [5.6402, 0.0, 0.0],
        [0.0, 5.6402, 0.0],
        [0.0, 0.0, 5.6402]
      ],
      'scaledPositions': [
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.5]
      ],
      'primitiveCell': [
        [0.0, 2.8201, 2.8201],
        [2.8201, 0.0, 2.8201],
        [2.8201, 2.8201, 0.0]
      ],
      'pbc': [true, true, true],
      'bonds': [[0, 1]]
    }
    viewer.current.load(bulk)
    viewer.current.resizeCanvasToHostElement()
    viewer.current.fitToCanvas()
    viewer.current.saveReset()
    viewer.current.reset()
  }, [])

  // These handle the viewer settings
  useEffect(() => {
    viewer.current.toggleBonds(showBonds)
  }, [showBonds])
  useEffect(() => {
    viewer.current.toggleLatticeParameters(showAxis)
  }, [showAxis])

  // Handles fullscreen mode. Needs to be done before render with
  // useLayoutEffect to avoid stuttering.
  useLayoutEffect(() => {
    if (viewer.current !== null) {
      viewer.current.resizeCanvasToHostElement()
      viewer.current.fitToCanvas()
      viewer.current.render()
    }
  }, [fullscreen])

  return (
    <Box className={classes.root}>
      <ClickAwayListener onClickAway={handleClickAway}>
        <Paper elevation={fullscreen ? 12 : 0} className={fullscreen ? [classes.container, classes.fullscreen].join(' ') : classes.container} ref={refContainer}>
          <div className={classes.header}>
            <IconButton className={classes.iconButton} onClick={handleOpenMenu}>
              <MoreVertIcon />
            </IconButton>
            <div className={classes.spacer}></div>
            <IconButton className={classes.iconButton} onClick={handleReset}>
              <ReplayIcon />
            </IconButton>
            <IconButton className={classes.iconButton} onClick={handleScreencapture}>
              <CameraAltIcon />
            </IconButton>
            <IconButton className={classes.iconButton} onClick={handleFullscreen}>
              <FullscreenIcon />
            </IconButton>
          </div>
          <Menu
            id='settings-menu'
            anchorEl={anchorEl}
            getContentAnchorEl={null}
            anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
            transformOrigin={{ vertical: 'top', horizontal: 'left' }}
            keepMounted
            open={open}
            onClose={handleCloseMenu}
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
                    checked={showAxis}
                    onChange={(event) => { setShowAxis(!showAxis) }}
                    color='primary'
                  />
                }
                label='Show axis'
              />
            </MenuItem>
          </Menu>
          <div className={classes.viewerCanvas} ref={refViewerCanvas}></div>
        </Paper>
      </ClickAwayListener>
    </Box>
  )
}

Structure.propTypes = {
  info: PropTypes.object,
  calcId: PropTypes.string
}
