import React, { useState, useEffect, useRef, useLayoutEffect } from 'react';
import PropTypes from 'prop-types'

import { makeStyles } from '@material-ui/core/styles';
import { ButtonGroup, Button, Card, Paper } from '@material-ui/core'
import IconButton from '@material-ui/core/IconButton';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import FullscreenIcon from '@material-ui/icons/Fullscreen';
import CameraAltIcon from '@material-ui/icons/CameraAlt';
import ReplayIcon from '@material-ui/icons/Replay';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Checkbox from '@material-ui/core/Checkbox';
import CheckBoxIcon from '@material-ui/icons/CheckBox';
import Box from '@material-ui/core/Box';

export function StructureViewer(props) {
  const [anchorEl, setAnchorEl] = React.useState(null);
  const open = Boolean(anchorEl);
  let viewerCanvas = useRef(null);

  let viewerOptions = {
    view: {
      autoResize: false,       // Automatically fit canvas to the root visualization element
      autoFit: false,          // Automatically set the zoom to fit the structure into the canvas
      fitMargin: 0.5,
    },
    controls: {
      enableZoom: true,       // Enable zoom with mouse wheel/pinch
      enablePan: true,        // Enable pan with mouse/touch
      enableRotate: true,     // Enable rotation
      zoomLevel: 1,           // Default zoom level
    },
    structure: {
      showParam: true,        // Show lattice parameters
      showCopies: false,      // Show periodic copies at cell boundary
      showShadows: true,      // Enable shadows: requires a bit more from GPU
      allowRepeat: false,     // 
      showCell: true,         // Show unit cell wireframe
      wrap: false,            // Wrap atoms to unit cell
      createLegend: false,    // Show atom labels
      showLegend: true,       // Show atom labels
      showOptions: true,      // Show the options toolbar
      showBonds: true,        // Show or hide bonds
      radiusScale: 0.8,       // Scaling factor for atomic radii
      bondScale: 1.2,         // Scaling factor for the automatic bond detection
      viewCenter: "COC",      // The rotation and view center position. Valid values are: "COP" (center of position), "COC" center of cell, or a custom position given as array of three cartesian coordinates.
    },
    style: {
      backgroundColor: [0xffffff, 1]
    }
  };

  // State
  let [options, setOptions] = useState(viewerOptions);
  let [viewer, setViewer] = useState(null);
  let [fullscreen, setFullscreen] = useState(false);

  // Default styles
  const useStyles = makeStyles({
    root: {
      width: "100%",
    },
    container: {
      boxSizing: "border-box",  // This should be a global default?
      padding: "0.5rem",
      display: "flex",
      height: "100%",
      width: "100%",
      flexDirection: "column",
      backgroundColor: "white",
    },
    fullscreen: {
      position: "fixed",
      zIndex: 2000,
      left: "50%",
      top: "50%",
      width: "80%",
      height: "50%",
      transform: "translate(-50%, -50%)",
    },
    header: {
      display: "flex",
      flexDirection: "row",
    },
    spacer: {
      flex: 1,
    },
    viewerCanvas: {
      flex: 1,
      minHeight: 0,  // added min-height: 0 to allow the item to shrink to fit inside the container.
    },
  });

  // Settings menu structure
  let showBonds = true;
  const optionMenu = [
    {label: "Show bonds", value: options.structure.showBonds, func: () => {
      options.structure.showBonds = !options.structure.showBonds;
      setOptions(options);
    }}
  ];

  const handleOpenMenu = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleCloseMenu = () => {
    setAnchorEl(null);
  };

  const handleFullscreen = () => {
    setFullscreen(!fullscreen);
  };

  const handleScreencapture = () => {
    viewer.takeScreenShot(props.calcId);
  };

  const handleReset = () => {
    viewer.reset();
    viewer.fitToCanvas();
    viewer.render();
  };

  // Run only on first render to initialize the viewer
  useEffect(() => {
    console.log("INITIALIZE VIEW");
    viewer = new window.matviewer.StructureViewer(viewerCanvas.current, viewerOptions);
    setViewer(viewer);
    var bulk = {
        "atomicNumbers": [11, 17, 11, 17, 11, 17, 11, 17],
        "cell": [
            [5.6402, 0.0, 0.0],
            [0.0, 5.6402, 0.0],
            [0.0, 0.0, 5.6402]
        ],
        "scaledPositions": [
            [0.0, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.0, 0.5]
        ],
        "primitiveCell": [
            [0.0, 2.8201, 2.8201],
            [2.8201, 0.0, 2.8201],
            [2.8201, 2.8201, 0.0]
        ],
        "pbc": [true, true, true],
        "bonds": [[0, 1]]
    };
    viewer.load(bulk);
    viewer.resizeCanvasToHostElement();
    viewer.fitToCanvas();
    viewer.saveReset();
    viewer.reset();
  }, []);

  // These handle the viewer settings
  useEffect(() => {
    viewer.toggleBonds(options.structure.showBonds);
  }, [options.structure.showBonds]);

  // Handles fullscreen mode
  useLayoutEffect(() => {
    if (viewer !== null) {
      viewer.resizeCanvasToHostElement();
      viewer.fitToCanvas();
      viewer.render();
    }
  }, [fullscreen]);

  const classes = useStyles(props);
  let containerClasses = fullscreen ? [classes.container, classes.fullscreen].join(" ") : classes.container;
  let elevation = fullscreen ? 12 : 0;

  return (
    <Box className={classes.root}>
      <Paper elevation={elevation} className={containerClasses}>
        <div className={classes.header}>
          <IconButton onClick={handleOpenMenu}>
            <MoreVertIcon />
          </IconButton>
          <div className={classes.spacer}></div>
          <IconButton onClick={handleReset}>
            <ReplayIcon />
          </IconButton>
          <IconButton onClick={handleScreencapture}>
            <CameraAltIcon />
          </IconButton>
          <IconButton onClick={handleFullscreen}>
            <FullscreenIcon />
          </IconButton>
        </div>
        <Menu
          id="settings-menu"
          anchorEl={anchorEl}
          getContentAnchorEl={null}
          anchorOrigin={{ vertical: "bottom", horizontal: "left" }}
          transformOrigin={{ vertical: "top", horizontal: "left" }}
          keepMounted
          open={open}
          onClose={handleCloseMenu}
        >
          {optionMenu.map((option) => (
            <MenuItem>
              <FormControlLabel
                control={
                  <Checkbox
                    checked={option.value}
                    onChange={option.func}
                    name="checkedB"
                    color="primary"
                  />
                }
                label={option.label}
              />
            </MenuItem>
          ))}
        </Menu>
        <div className={classes.viewerCanvas} ref={viewerCanvas}></div>
      </Paper>
    </Box>
  )
}

StructureViewer.propTypes = {
  info: PropTypes.object
}
