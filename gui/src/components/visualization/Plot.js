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
import { useHistory } from 'react-router-dom'

import {
  Typography,
  Box
} from '@material-ui/core'
import {
  Fullscreen,
  FullscreenExit,
  CameraAlt,
  Replay,
  ViewList
} from '@material-ui/icons'
import Floatable from './Floatable'
import Placeholder from '../visualization/Placeholder'
import Actions from '../Actions'
import Plotly from 'plotly.js-cartesian-dist-min'
import clsx from 'clsx'
import { mergeObjects } from '../../utils'

export default function Plot({
  data,
  layout,
  config,
  floatTitle,
  metaInfoLink,
  capture,
  aspectRatio,
  fixedMargins,
  className,
  classes,
  onRelayout,
  onRedraw,
  onRelayouting,
  onHover,
  onReset,
  layoutSubject
}) {
  // States
  const [float, setFloat] = useState(false)
  const [captureSettings, setCaptureSettings] = useState()
  const firstRender = useRef(true)
  const [margins, setMargins] = useState()
  const [loading, setLoading] = useState(true)
  const history = useHistory()

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
  const styles = useStyles({classes: classes})

  // Set the final layout
  const finalLayout = useMemo(() => {
    let defaultLayout = {
      dragmode: 'pan',
      hovermode: false,
      showlegend: false,
      autosize: true,
      // There is extra space reserved for the top and bottom margins so that
      // automargin does not make the plot jump around too much.
      margin: {
        l: 30,
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
        automargin: false,
        autorange: true,
        linecolor: '#333',
        linewidth: 1,
        mirror: true,
        ticks: 'outside',
        showline: true,
        fixedrange: true,
        title: {
          standoff: 10,
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
        autorange: true,
        linecolor: '#333',
        linewidth: 1,
        mirror: true,
        ticks: 'outside',
        showline: true,
        title: {
          standoff: 10,
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
    const finalLayoutResult = mergeObjects(layout, defaultLayout)
    // Override automargin settings
    if (fixedMargins && margins) {
      finalLayoutResult.yaxis.automargin = false
      finalLayoutResult.xaxis.automargin = false
      finalLayoutResult.margin.l = margins.l
      finalLayoutResult.margin.t = margins.t
    }
    return finalLayoutResult
  }, [fixedMargins, layout, margins])

  // Save the initial layout as reset
  const finalResetLayout = useMemo(() => {
    const resLayout = cloneDeep(finalLayout)
    return resLayout
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [margins])

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
    if (!firstRender.current) {
      // console.log('Plotting new data')
      Plotly.react(canvasRef.current, data, finalLayout, finalConfig)
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [firstRender, canvasRef, data, finalConfig])

  // Relayout plots when layout updated. useLayoutEffect is used so that React
  // rendering and Plotly rendering are synced.
  useLayoutEffect(() => {
    if (!firstRender.current && canvasRef?.current && finalLayout) {
      // console.log('Relayouting')
      // console.log(finalLayout)
      Plotly.relayout(canvasRef.current, finalLayout)
    }
  }, [firstRender, canvasRef, finalLayout])

  // Only called once after the first render. Sets the reference for determining
  // if an initial render has been done. Using a reference instead of a state
  // avoids a rerender.
  useLayoutEffect(() => {
    if (firstRender.current) {
      firstRender.current = false
    }
  }, [])

  // Captures the first margin value from the plot.
  useLayoutEffect(() => {
    if (canvasRef.current && fixedMargins) {
      // console.log('Figuring out new margins')
      try {
        // Get the element which explicitly stores the computed margin, and save
        // these values on the layout directly.
        const group = canvasRef.current.children[0].children[0].children[0].children[2].children[0].children[0]
        const l = parseInt(group.getAttribute('x'))
        const t = parseInt(group.getAttribute('y'))
        setMargins({l: l, t: t})
      } catch (e) {
        console.log('Could not determine the margin values.')
      }
    }
  }, [canvasRef, data, fixedMargins])

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
  if (metaInfoLink) {
    actions.push({tooltip: 'View data in the archive', onClick: () => { history.push(metaInfoLink) }, content: <ViewList/>})
  }

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
  fixedMargins: PropTypes.bool, // Whether to automatically update margins beyond the first render.
  className: PropTypes.string,
  classes: PropTypes.string,
  floatTitle: PropTypes.string, // The title of the plot shown in floating mode
  metaInfoLink: PropTypes.string, // A link to the metainfo definition
  onRelayout: PropTypes.func,
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
  floatTitle: '',
  fixedMargins: true
}
