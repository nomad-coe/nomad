import React, { useState, useEffect, useMemo, useRef, useCallback } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import { cloneDeep } from 'lodash'

import {
  IconButton,
  Tooltip,
  Typography
} from '@material-ui/core'
import {
  MoreVert,
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay
} from '@material-ui/icons'
import Floatable from './Floatable'
import Plotly from 'plotly.js-cartesian-dist-min'
import clsx from 'clsx'
import { mergeObjects } from '../../utils'

export default function Plot({data, layout, config, menu, floatTitle, capture, aspectRatio, className, classes, onRelayout, onAfterPlot, onRedraw, onRelayouting}) {
  // States
  const [float, setFloat] = useState(false)
  const [initialLayout, setInitialLayout] = useState()
  const [captureSettings, setCaptureSettings] = useState()
  const firstUpdate = useRef(true)
  const dataInitialized = useRef(false)

  React.useEffect(() => {
    let defaultCapture = {
      format: 'png',
      width: 1024,
      height: 960 / aspectRatio,
      filename: 'plot'
    }
    let settings = mergeObjects(capture, defaultCapture)
    setCaptureSettings(settings)
  }, [capture, aspectRatio])

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      header: {
        paddingRight: theme.spacing(1),
        display: 'flex',
        flexDirection: 'row',
        zIndex: 1
      },
      root: {
      },
      spacer: {
        flex: 1
      },
      iconButton: {
        backgroundColor: 'white'
      }
    }
  })

  const style = useStyles(classes)

  // Set the final menu
  const finalMenu = useMemo(() => {
    let defaultMenu = {
      reset: {
        visible: true,
        disabled: false
      },
      fullscreen: {
        visible: true,
        disabled: false
      },
      capture: {
        visible: true,
        disabled: false
      },
      dropdown: {
        visible: false,
        disabled: false,
        items: undefined
      }
    }
    return mergeObjects(menu, defaultMenu)
  }, [menu])

  // Set the final layout
  const finalLayout = useMemo(() => {
    let defaultLayout = {
      dragmode: 'pan',
      hovermode: false,
      showlegend: false,
      autosize: true,
      margin: {
        l: 60,
        r: 20,
        t: 20,
        b: 50
      },
      xaxis: {
        linecolor: '#333',
        linewidth: 1,
        mirror: true,
        ticks: 'outside',
        showline: true,
        fixedrange: true,
        title: {
          font: {
            family: 'Titillium Web,sans-serif',
            size: 16,
            color: '#333'
          },
          tickfont: {
            family: 'Titillium Web,sans-serif',
            size: 14,
            color: '#333'
          }
        }
      },
      yaxis: {
        automargin: true,
        linecolor: '#333',
        linewidth: 1,
        mirror: true,
        ticks: 'outside',
        showline: true,
        title: {
          font: {
            family: 'Titillium Web,sans-serif',
            size: 16,
            color: '#333'
          }
        },
        tickfont: {
          family: 'Titillium Web,sans-serif',
          size: 14,
          color: '#333'
        }
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout])

  // Set the final config
  const finalConfig = useMemo(() => {
    let defaultConfig = {
      scrollZoom: true,
      displayModeBar: false
    }
    return mergeObjects(config, defaultConfig)
  }, [config])

  // Initialize the plot object on first render
  useEffect(() => {
    Plotly.newPlot(canvasRef.current, data, finalLayout, finalConfig)
    if (firstUpdate.current) {
      firstUpdate.current = false
    }
    if (data) {
      setInitialLayout(cloneDeep(finalLayout))
      dataInitialized.current = true
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // In order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const canvasRef = useCallback(node => {
    canvasRef.current = node
    if (node === null) {
      return
    }
    Plotly.react(canvasRef.current, data, finalLayout, finalConfig)
  }, [data, finalLayout, finalConfig])

  // Update data
  useEffect(() => {
    if (firstUpdate.current) {
      return
    }
    Plotly.react(canvasRef.current, data, finalLayout, finalConfig)
    if (data) {
      if (!dataInitialized.current) {
        setInitialLayout(cloneDeep(finalLayout))
        dataInitialized.current = true
      }
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data, finalConfig])

  // Update layout
  useEffect(() => {
    if (firstUpdate.current) {
      return
    }
    Plotly.relayout(canvasRef.current, finalLayout)
  }, [canvasRef, finalLayout])

  // For resetting the view. We need to pass a deep clone of the initialLayout,
  // as Plotly will modify the given object in-place
  const handleReset = useCallback(() => {
    if (initialLayout) {
      Plotly.relayout(canvasRef.current, cloneDeep(initialLayout))
    }
  }, [canvasRef, initialLayout])

  // Handles plot capturing
  const handleCapture = useCallback(() => {
    Plotly.downloadImage(canvasRef.current, captureSettings)
  }, [canvasRef, captureSettings])

  return (
    <Floatable className={clsx(className, style.root)} float={float} onFloat={() => setFloat(!float)} aspectRatio={aspectRatio}>
      <div className={style.header}>
        {float && <Typography variant="h6">{floatTitle}</Typography>}
        <div className={style.spacer}></div>
        { finalMenu.reset.visible === true
          ? <Tooltip title="Reset view">
            <IconButton className={style.iconButton} onClick={handleReset} disabled={finalMenu.reset.disabled}> <Replay />
            </IconButton>
          </Tooltip>
          : ''
        }
        { finalMenu.fullscreen.visible === true
          ? <Tooltip
            title="Toggle fullscreen">
            <IconButton className={style.iconButton} onClick={() => setFloat(!float)} disabled={finalMenu.fullscreen.disabled}>
              {float ? <FullscreenExit /> : <Fullscreen />}
            </IconButton>
          </Tooltip>
          : ''
        }
        { finalMenu.capture.visible === true
          ? <Tooltip title="Capture image">
            <IconButton className={style.iconButton} onClick={handleCapture} disabled={finalMenu.capture.disabled}>
              <CameraAlt />
            </IconButton>
          </Tooltip>
          : ''
        }
        { finalMenu.dropdown.visible === true
          ? <Tooltip title="Options">
            <IconButton className={style.iconButton} onClick={() => {}} disabled={finalMenu.dropdown.disabled}>
              <MoreVert />
            </IconButton>
          </Tooltip>
          : ''
        }
      </div>
      <div ref={canvasRef} style={{width: '100%', height: '100%'}}></div>
    </Floatable>
  )
}

Plot.propTypes = {
  data: PropTypes.array, // Plotly.js data object
  layout: PropTypes.object, // Plotly.js layout object
  config: PropTypes.object, // Plotly.js config object
  menu: PropTypes.object, // Menu settings
  capture: PropTypes.object, // Capture settings
  aspectRatio: PropTypes.number, // Fixed aspect ratio for the viewer canvas
  className: PropTypes.string,
  classes: PropTypes.string,
  floatTitle: PropTypes.string, // The title of the plot shown in floating mode
  onRelayout: PropTypes.func,
  onAfterPlot: PropTypes.func,
  onRedraw: PropTypes.func,
  onRelayouting: PropTypes.func
}
Plot.defaultProps = {
  aspectRatio: 9 / 16,
  floatTitle: ''
}
