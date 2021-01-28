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
import React, { useState, useEffect, useLayoutEffect, useMemo, useRef, useCallback } from 'react'
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

export default function Plot({data, layout, config, menu, floatTitle, capture, aspectRatio, className, classes, onRelayout, onAfterPlot, onRedraw, onRelayouting, onHover, onReset}) {
  // States
  const [float, setFloat] = useState(false)
  const [captureSettings, setCaptureSettings] = useState()
  const firstUpdate = useRef(true)

  useEffect(() => {
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
      title: {
        font: {
          family: 'Titillium Web,sans-serif'
        }
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

  // Save the initial layout as reset
  const finalResetLayout = useMemo(() => {
    return cloneDeep(finalLayout)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

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
      // Attach events when the canvas is first created
      if (onRelayouting) {
        canvasRef.current.on('plotly_relayouting', onRelayouting)
      }
      if (onRedraw) {
        canvasRef.current.on('plotly_redraw', onRedraw)
      }
      if (onRelayout) {
        canvasRef.current.on('plotly_relayout', onRelayout)
      }
      if (onHover) {
        canvasRef.current.on('plotly_hover', onHover)
      }
      firstUpdate.current = false
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // This callback redraws the plot whenever the canvas element changes. In
  // order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const canvasRef = useCallback(node => {
    // Do nothing if the canvas has not actually changed.
    if (node === null || canvasRef.current === node) {
      return
    }

    // Redraw if moving to a new canvas. Also re-attach events.
    if (canvasRef.current !== undefined) {
      let oldLayout = canvasRef.current.layout
      let oldData = canvasRef.current.data
      // Update canvas and redraw on it
      Plotly.react(node, oldData, oldLayout, finalConfig)

      // Re-attach events whenever the canvas changes
      if (onRelayouting) {
        node.on('plotly_relayouting', onRelayouting)
      }
      if (onRedraw) {
        node.on('plotly_redraw', onRedraw)
      }
      if (onRelayout) {
        node.on('plotly_relayout', onRelayout)
      }
      if (onHover) {
        node.on('plotly_hover', onHover)
      }
    }

    // Update canvas element
    canvasRef.current = node

  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Update plots when data or config is updated. useLayoutEffect is used so
  // that React rendering and Plotly rendering are synced.
  useLayoutEffect(() => {
    if (!firstUpdate.current) {
      Plotly.react(canvasRef.current, data, finalLayout, finalConfig)
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [canvasRef, data, finalConfig])

  // Relayout plots when layout updated. useLayoutEffect is used so that React
  // rendering and Plotly rendering are synced.
  useLayoutEffect(() => {
    if (!firstUpdate.current && canvasRef?.current && finalLayout) {
      Plotly.relayout(canvasRef.current, finalLayout)
    }
  }, [canvasRef, finalLayout])

  // For resetting the view.
  const handleReset = useCallback(() => {
    if (canvasRef?.current && finalResetLayout) {
      Plotly.relayout(canvasRef.current, mergeObjects(finalResetLayout, canvasRef.current.layout))
      if (onReset) {
        onReset()
      }
    }
  }, [canvasRef, finalResetLayout, onReset])

  // Handles plot capturing
  const handleCapture = useCallback(() => {
    Plotly.downloadImage(canvasRef.current, captureSettings)
  }, [canvasRef, captureSettings])

  return (
    <Floatable className={clsx(className, style.root)} float={float} onFloat={() => setFloat(!float)} aspectRatio={aspectRatio}>
      {float && <Typography variant="h6">{floatTitle}</Typography>}
      <div ref={canvasRef} style={{width: '100%', height: '100%'}}></div>
      <div className={style.header}>
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
  onRelayouting: PropTypes.func,
  onHover: PropTypes.func,
  onReset: PropTypes.func
}
Plot.defaultProps = {
  aspectRatio: 9 / 16,
  floatTitle: ''
}
