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
  Typography,
  Box
} from '@material-ui/core'
import {
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay
} from '@material-ui/icons'
import Floatable from './Floatable'
import Placeholder from '../visualization/Placeholder'
import Actions from '../Actions'
import Plotly from 'plotly.js-cartesian-dist-min'
import clsx from 'clsx'
import { mergeObjects } from '../../utils'

export default function Plot({data, layout, config, floatTitle, capture, aspectRatio, className, classes, onRelayout, onAfterPlot, onRedraw, onRelayouting, onHover, onReset, layoutSubject}) {
  // States
  const [float, setFloat] = useState(false)
  const [captureSettings, setCaptureSettings] = useState()
  const firstUpdate = useRef(true)
  const [loading, setLoading] = useState(true)

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
      placeHolder: {
        left: 0,
        right: 0,
        position: 'absolute'
      },
      floatable: {
        visibility: loading ? 'hidden' : 'visible'
      },
      spacer: {
        flex: 1
      },
      iconButton: {
        backgroundColor: 'white',
        marginLeft: theme.spacing(1)
      }
    }
  })

  const styles = useStyles(classes)

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
        autorange: true,
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
      displayModeBar: false,
      showTips: false
    }
    return mergeObjects(config, defaultConfig)
  }, [config])

  // This callback redraws the plot whenever the canvas element changes. In
  // order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const canvasRef = useCallback(node => {
    // Do nothing if the canvas has not actually changed or data is not ready
    if (node === null || canvasRef.current === node || !data) {
      return
    }

    // When the canvas reference is instantiated for the first time, create a
    // new plot.
    if (canvasRef.current === undefined) {
      Plotly.newPlot(node, data, finalLayout, finalConfig)
      if (firstUpdate.current) {
        firstUpdate.current = false
      }
      setLoading(false)
    // When the reference changes for the second time, react instead to save
    // some time
    } else {
      let oldLayout = canvasRef.current.layout
      let oldData = canvasRef.current.data
      Plotly.react(node, oldData, oldLayout, finalConfig)
    }

    // (Re-)attach events whenever the canvas changes
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

    // Subscribe to the layout change publisher if one is given
    if (layoutSubject) {
      layoutSubject.subscribe(layout => {
        let oldLayout = canvasRef.current.layout
        Plotly.relayout(canvasRef.current, mergeObjects(layout, oldLayout))
      })
    }

    // Update canvas element
    canvasRef.current = node
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [data, onHover, onRedraw, onRelayout, onRelayouting])

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

  // List of actionable buttons for the plot
  const actions = [
    {tooltip: 'Reset view', onClick: handleReset, content: <Replay/>},
    {tooltip: 'Toggle fullscreen', onClick: () => setFloat(!float), content: float ? <FullscreenExit/> : <Fullscreen/>},
    {tooltip: 'Capture image', onClick: handleCapture, content: <CameraAlt/>}
  ]

  // Even if the plots are still loading, all the html elements need to be
  // placed in the DOM. During loading, they are placed underneath the
  // palceholder with visibility=hidden. This way Plotly still has access to
  // these HTML nodes and their sizes when the plots are loading.
  return <Box className={clsx(className, styles.root)} position='relative' width='100%'>
    {loading && <Placeholder className={styles.placeHolder} variant="rect" aspectRatio={aspectRatio}></Placeholder>}
    <Floatable className={styles.floatable} float={float} onFloat={() => setFloat(!float)} aspectRatio={aspectRatio}>
      {float && <Typography variant="h6">{floatTitle}</Typography>}
      <div ref={canvasRef} style={{width: '100%', height: '100%'}}></div>
      <div className={styles.header}>
        <Actions actions={actions}></Actions>
      </div>
    </Floatable>
  </Box>
}

Plot.propTypes = {
  data: PropTypes.array, // Plotly.js data object
  layout: PropTypes.object, // Plotly.js layout object
  config: PropTypes.object, // Plotly.js config object
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
  onReset: PropTypes.func,
  /**
   * A RxJS Subject for efficient, non-persistent, layout changes that bypass
   * rendering of the component. Should send messages that contain the new
   * layout object.
  */
  layoutSubject: PropTypes.any
}
Plot.defaultProps = {
  aspectRatio: 9 / 16,
  floatTitle: ''
}
