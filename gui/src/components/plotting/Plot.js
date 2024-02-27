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
import React, {useState, useLayoutEffect, useMemo, useRef, useCallback, forwardRef, useImperativeHandle, useEffect} from 'react'
import PropTypes from 'prop-types'
import { useResizeDetector } from 'react-resize-detector'
import { isNil, isEmpty, cloneDeep } from 'lodash'
import { useHistory } from 'react-router-dom'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Typography,
  Box,
  Tooltip
} from '@material-ui/core'
import {
  CameraAlt,
  Fullscreen,
  FullscreenExit,
  Replay,
  ViewList,
  Warning
} from '@material-ui/icons'
import Floatable from '../visualization/Floatable'
import Placeholder from '../visualization/Placeholder'
import NoData from '../visualization/NoData'
import { Actions, Action } from '../Actions'
import Plotly from 'plotly.js-dist-min'
import clsx from 'clsx'
import { mergeObjects } from '../../utils'

/**
 * Component that produces a highly customized Plotly.js graph. We are using the
 * vanilla JS Plotly library to achieve the most customizability.
 */
const useStyles = makeStyles((theme) => {
  return {
    root: {
      width: '100%',
      height: '100%',
      position: 'relative'
    },
    footer: {
      paddingRight: theme.spacing(1),
      display: 'flex',
      flexDirection: 'row',
      zIndex: 1
    },
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
    canvasContainer: {
      flex: 1,
      minHeight: 0 // added min-height: 0 to allow the item to shrink to fit inside the container.
    },
    placeholderRoot: {
      left: 0,
      right: 0,
      position: 'absolute'
    },
    floatable: {
      width: '100%',
      height: '100%'
    },
    spacer: {
      flex: 1
    },
    warning: {
      color: theme.palette.warning.main,
      position: 'absolute',
      zIndex: 1,
      left: '2px',
      top: '20px'
    }
  }
})
const Plot = React.memo(forwardRef(({
  data,
  layout,
  config,
  floatTitle,
  metaInfoLink,
  capture,
  fixedMargins,
  className,
  classes,
  placeHolderStyle,
  noDataStyle,
  warning,
  onRelayout,
  onRedraw,
  onRelayouting,
  onHover,
  onReset,
  onSelected,
  onSelecting,
  onDeselect,
  onClick,
  onWebglContextLost,
  layoutSubject,
  disableDefaultActions,
  actions,
  throttleResize,
  'data-testid': testID
}, canvas) => {
  const [float, setFloat] = useState(false)
  const theme = useTheme()
  const firstRender = useRef(true)
  const [ratio, setRatio] = useState(1)
  const [canvasNode, setCanvasNode] = useState()
  const attach = useRef()
  const [sizeReady, setSizeReady] = useState(false)
  const [loading, setLoading] = useState(true)
  const margins = useRef()
  const history = useHistory()
  const { height, width, ref } = useResizeDetector()
  const canvasSize = useRef({})
  const layoutRecoverRef = useRef()
  const resetInfo = useRef()

  // The image capture settings
  const captureSettings = useMemo(() => {
    const maxSize = 1280
    let width, height
    if (ratio > 1) {
      width = maxSize
      height = maxSize / ratio
    } else {
      width = maxSize * ratio
      height = maxSize
    }
    const defaultCapture = {
      format: 'png',
      width: width,
      height: height,
      filename: 'plot'
    }
    const settings = mergeObjects(capture, defaultCapture)
    return settings
  }, [capture, ratio])

  const styles = useStyles({classes: classes})
  if (noDataStyle) {
    styles.nodata = noDataStyle
  }
  if (placeHolderStyle) {
    styles.placeholder = placeHolderStyle
  }

  // This callback redraws the plot whenever the canvas element changes. In
  // order to properly detect changes in a reference, a reference callback is
  // used. This is the recommended way to monitor reference changes as a simple
  // useRef is not guaranteed to update:
  // https://reactjs.org/docs/hooks-faq.html#how-can-i-measure-a-dom-node
  const canvasRef = useCallback(node => {
    if (node) {
      if (!firstRender.current) {
        attach.current = true
        // We use a dummy state to request a redraw on the new canvas
        setCanvasNode(node)
      }
      // We add a NATIVE wheel event handler that prevents the wheel events sent
      // by Plotly from bubbling up. Notice that the SyntheticEvents added by
      // React cannot be used, since they do not play well with the original
      // native events created by Plotly (in particular: SyntheticEvents are all
      // passive, so calling preventDefault on them will cause problems on
      // browsers that support passive events.)
      node.addEventListener('wheel', (e) => {
        e.preventDefault()
        e.stopPropagation()
      })
    }
    canvasRef.current = node
  }, [])

  // Set the final layout. It is a combination of a default layout, the layout
  // set by the user and some properties of the curretly used layout.
  const finalLayout = useMemo(() => {
    const withTitle = layout?.title?.text || layout?.template?.title?.text || layout?.annotations?.some(item => item?.text)
    const defaultLayout = {
      dragmode: 'pan',
      hovermode: false,
      showlegend: false,
      autosize: false,
      // There is extra space reserved for the top and bottom margins so that
      // automargin does not make the plot jump around too much.
      margin: {
        l: theme.spacing(4),
        r: theme.spacing(1.5),
        t: theme.spacing(withTitle ? 5 : 1),
        b: theme.spacing(6)
      },
      title: {
        font: {
          family: theme.typography.fontFamily
        },
        yanchor: 'middle'
      },
      legend: {
        bgcolor: 'rgba(255, 255, 255, 0.9)',
        font: {
          family: theme.typography.fontFamily,
          size: 14
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
            family: theme.typography.fontFamily,
            size: 16,
            color: '#333'
          },
          tickfont: {
            family: theme.typography.fontFamily,
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
            family: theme.typography.fontFamily,
            size: 16,
            color: '#333'
          }
        },
        tickfont: {
          family: 'Titillium Web,sans-serif',
          size: 14,
          color: '#333'
        }
      },
      yaxis2: {
        automargin: true,
        autorange: true,
        linecolor: '#333',
        linewidth: 1,
        mirror: true,
        ticks: 'outside',
        showline: true,
        title: {
          standoff: 8,
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
    const newLayout = mergeObjects(layout, defaultLayout)
    return newLayout
  }, [layout, theme])

  // Set the final config
  const finalConfig = useMemo(() => {
    const defaultConfig = {
      scrollZoom: true,
      displayModeBar: false,
      showTips: false,
      responsive: true
    }
    return mergeObjects(config, defaultConfig)
  }, [config])

  // Saves reset options according to the given layout.
  useEffect(() => {
    resetInfo.current = {
        rangeX: finalLayout?.xaxis?.range,
        rangeY: finalLayout?.yaxis?.range,
        rangeY2: finalLayout?.yaxis2?.range,
        margin: finalLayout?.margin
    }
  }, [finalLayout])

  // Used to get the currently active layout. Because of Floatable, the plot
  // layout is lost upon the plot being put on another div. This is why we need
  // a backup that is stored when the final plot config is ready. Notice that we
  // have to keep parts of the old layout, as it stores e.g. range etc. that
  // need to be persistent.
  const getLayout = useCallback((update) => {
      const oldLayout = cloneDeep(canvasRef.current.layout || layoutRecoverRef.current || {})
      if (!oldLayout.xaxis) oldLayout.xaxis = {}
      if (!oldLayout.yaxis) oldLayout.yaxis = {}
      if (!oldLayout.yaxis2) oldLayout.yaxis2 = {}
      const newLayout = mergeObjects(update, oldLayout)
      // When fixed margins are set, the margins should be automatically set and
      // then fixed in a secondary update.
      if (fixedMargins) {
        newLayout.xaxis.automargin = true
        newLayout.yaxis.automargin = true
        newLayout.yaxis2.automargin = true
        newLayout.margin = resetInfo.current.margin
      }
      // The current size overrides any size that is set
      newLayout.width = canvasSize.current?.width || newLayout.width
      newLayout.height = canvasSize.current?.height || newLayout.height
      // The range from updated layout should be used if it's set, otherwise use
      // old range.
      if (isEmpty(newLayout.xaxis.range)) newLayout.xaxis.range = oldLayout.xaxis.range
      if (isEmpty(newLayout.yaxis.range)) newLayout.yaxis.range = oldLayout.yaxis.range
      if (isEmpty(newLayout.yaxis2.range)) newLayout.yaxis2.range = oldLayout.yaxis2.range
      return newLayout
  }, [canvasRef, fixedMargins])

  // If fixedMargins = true, captures the automatically set margin value from
  // the plot and performs a layout update that fixes them.
  const fixMargins = useCallback(() => {
    if (fixedMargins && canvasRef.current) {
      try {
        // Get the element which explicitly stores the computed margin and
        // perform a relayout with these values.
        const currentLayout = cloneDeep(canvasRef.current.layout)
        currentLayout.xaxis.automargin = false
        currentLayout.yaxis.automargin = false
        currentLayout.yaxis2.automargin = false
        if (!currentLayout.margin) currentLayout.margin = {}
        const subPlots = canvasRef.current.children[0].children[0].children[0].children[2].children
        let [xMin, xMax] = [undefined, undefined]
        for (const subPlot of subPlots) {
          const canvas = subPlot.children[0]
          const left = parseInt(canvas.getAttribute('x'))
          const right = left + parseInt(canvas.getAttribute('width'))
          xMin = !xMin || left < xMin ? left : xMin
          xMax = !xMax || right > xMax ? right : xMax
        }
        const totalWidth = canvasRef.current.clientWidth
        currentLayout.margin.l = xMin
        currentLayout.margin.r = totalWidth - xMax
        margins.current = cloneDeep(currentLayout.margin)
        Plotly.relayout(canvasRef.current, currentLayout)
      } catch (e) {
        console.log('Could not determine the margin values.', e)
      }
    }
  }, [canvasRef, fixedMargins])

  // Handle canvas resize
  useLayoutEffect(() => {
    if (width && height) {
      setRatio(width / height)
      canvasSize.current = {width, height}
      if (firstRender.current) {
        setSizeReady(true)
      } else {
        // The updates are throttled by using requestAnimationFrame: there is no
        // sense in trying to update beyond what the browser can render TODO:
        // Ideally we would not need to call fixMargins here, as a simple resize
        // should not trigger any change in the margins. The call is curently
        // necesary due to the way the plot updates when going fullscreen:
        // without fixmargins the view is jumbled.
        const layout = getLayout()
        if (throttleResize) {
          window.requestAnimationFrame(() => {
            // Notice that fixMargins has to be done within this callback
            // together with the relayout
            Plotly.relayout(canvasRef.current, layout)
            fixMargins()
          })
        } else {
          Plotly.relayout(canvasRef.current, layout)
          fixMargins()
        }
      }
    }
  }, [canvasRef, width, height, throttleResize, fixMargins, getLayout])

  // Update plots when data, config or rendering canvas is updated.
  // useLayoutEffect is used so that React rendering and Plotly rendering are
  // synced.
  useLayoutEffect(() => {
    if (!data || !sizeReady) {
      if (data === false) setLoading(false)
      return
    }
    // The plot is created on the first render
    const layout = getLayout(finalLayout)
    if (firstRender.current) {
      Plotly.newPlot(canvasRef.current, data, layout, finalConfig)
      attach.current = true
      firstRender.current = false

      // Subscribe to the layout change publisher if one is given
      if (layoutSubject) {
        layoutSubject.subscribe(layout => {
          // The updates are throttled by using requestAnimationFrame: there is
          // no sense in trying to update beyond what the browser can render
          const newLayout = mergeObjects(layout, canvasRef.current.layout)
          window.requestAnimationFrame(() => {
            Plotly.relayout(canvasRef.current, newLayout)
          })
        })
      }
      setLoading(false)
    // On subsequent updates to data or layout, we call Plotly.react to update
    // the plot.
    } else {
      Plotly.react(canvasRef.current, data, layout, finalConfig)
    }

    fixMargins()

    // Upon first render or changing the plot DOM element the events are (re-)attached.
    if (attach.current === true) {
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
      if (onSelected) {
        canvasRef.current.on('plotly_selected', onSelected)
      }
      if (onSelecting) {
        canvasRef.current.on('plotly_selecting', onSelecting)
      }
      if (onDeselect) {
        canvasRef.current.on('plotly_deselect', onDeselect)
      }
      if (onClick) {
        canvasRef.current.on('plotly_click', onClick)
      }
      if (onWebglContextLost) {
        canvasRef.current.on('plotly_webglcontextlost', onWebglContextLost)
      }
      attach.current = false
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [layoutSubject, fixedMargins, canvasNode, firstRender, data, finalConfig, finalLayout, sizeReady, canvasRef, onRelayout, onRelayouting, fixMargins, getLayout])

  // If the callbacks change, we need to update them in the canvas as well.
  // TODO: Do this for other callbacks as well.
  useEffect(() => {
    if (onSelected && canvasRef?.current?.on) {
      canvasRef.current.on('plotly_selected', onSelected)
    }
  }, [onSelected, canvasRef])

  // For resetting the view.
  const handleReset = useCallback(() => {
    if (canvasRef?.current) {
      // If a fixed range was set in the original layout, use it. Otherwise
      // autorange the view.
      const hasX = !isNil(resetInfo.current.rangeX)
      const hasY = !isNil(resetInfo.current.rangeY)
      const hasY2 = !isNil(resetInfo.current.rangeY2)
      const oldLayout = getLayout()
      const layout = {
        ...oldLayout,
        xaxis: {...oldLayout.xaxis, autorange: !hasX, range: hasX ? [...resetInfo.current.rangeX] : oldLayout.xaxis.range},
        yaxis: {...oldLayout.yaxis, autorange: !hasY, range: hasY ? [...resetInfo.current.rangeY] : oldLayout.yaxis.range},
        yaxis2: {...oldLayout.yaxis2, autorange: !hasY2, range: hasY2 ? [...resetInfo.current.rangeY2] : oldLayout.yaxis2.range}
      }
      Plotly.relayout(canvasRef.current, layout)
      fixMargins()
      onReset && onReset()
    }
  }, [canvasRef, onReset, fixMargins, getLayout])

  // Saves the current layout for recovering the plot
  const saveLayout = useCallback(() => {
    layoutRecoverRef.current = cloneDeep(canvasRef?.current?.layout)
  }, [canvasRef])

  // Handles float change
  const handleFloat = useCallback(() => {
    saveLayout()
    setFloat(!float)
  }, [float, saveLayout])

  // Handles plot capturing
  const handleCapture = useCallback(() => {
    Plotly.downloadImage(canvasRef.current, captureSettings)
  }, [canvasRef, captureSettings])

  // Listen to keyboard events to deselect.
  useEffect(() => {
    function handleDeselect(event) {
      if (['Escape', 'Delete'].includes(event.key)) {
        const layout = canvasRef.current.layout
        delete layout.selections
        const data = canvasRef.current.data
        data.forEach(element => {
          delete element.selectedpoints
        })
        Plotly.react(canvasRef.current, data, layout)
        onDeselect && onDeselect()
      }
    }
    window.addEventListener('keydown', handleDeselect)
    return () => window.removeEventListener('keydown', handleDeselect)
  }, [canvasRef, onDeselect])

  useImperativeHandle(canvas, () => ({
    reset() {
      handleReset()
    },
    toggleFullscreen() {
      handleFloat()
    },
    capture() {
      handleCapture()
    },
    saveLayout() {
      saveLayout()
    },
    relayout(callback) {
      if (canvasRef.current?.layout) {
        const newLayout = callback(canvasRef.current.layout)
        Plotly.relayout(canvasRef.current, newLayout)
      }
    }
  }))

  // Decides the possible overlay to display
  const overlay = useMemo(() => {
    if (loading) {
      return <Placeholder
        margin={0}
        className={styles.placeholderRoot}
        classes={{skeleton: placeHolderStyle}}
        variant="rect"
        data-testid={`${testID}-placeholder`}
      />
    } else if (data === false) {
      return <NoData classes={{placeholder: styles.nodata}} data-testid={`${testID}-nodata`}/>
    }
    return null
  }, [loading, data, styles, placeHolderStyle, testID])

  // Determine the final list of actions to show
  const allActions = useMemo(() => {
    const defaultActions = disableDefaultActions
      ? null
      : <>
      <Action key="reset" tooltip='Reset view' onClick={handleReset}>
        <Replay/>
      </Action>
      <Action key="fullscreen" tooltip='Toggle fullscreen' onClick={handleFloat}>
        {float ? <FullscreenExit/> : <Fullscreen/>}
      </Action>
      <Action key="capture" tooltip='Capture image' onClick={handleCapture}>
        <CameraAlt/>
      </Action>
      {metaInfoLink && <Action key="archive" tooltip='View data in the archive' onClick={() => { history.push(metaInfoLink) }}>
        <ViewList/>
      </Action>}
    </>
    return (defaultActions || actions)
      ? <div className={styles.footer}>
          <Actions>{defaultActions}{actions}</Actions>
        </div>
      : null
  }, [actions, disableDefaultActions, float, handleCapture, handleFloat, handleReset, history, metaInfoLink, styles])

  // Even if the plots are still loading, all the html elements need to be
  // placed in the DOM. During loading, they are placed underneath the
  // placeholder with visibility=hidden. This way Plotly still has access to
  // these HTML nodes and their sizes when the plots are loading.
  return <Box
    className={clsx(styles.root, className)}
    data-testid={testID}
  >
    <Floatable
      className={styles.floatable}
      float={float}
      onFloat={handleFloat}
    >
      <div className={styles.container}>
        {/* Overlay that is displayed on top of the viewer */}
        <div className={styles.content}>
          {overlay}
        </div>
        {/* This contains the actual visualization. Needs to be always in the DOM so we can
        load the structure into it properly. */}
        <div
          className={styles.content}
          style={{visibility: overlay ? 'hidden' : 'visible'}}
        >
          {float && <Typography variant="h6">{floatTitle}</Typography>}
          <div className={styles.canvasContainer} ref={ref}>
            {/* Note that we need to apply "style" to the canvas instead of "className" */}
            <div ref={canvasRef} style={{width: '100%', height: '100%', position: 'relative'}}>
              {warning && <Tooltip title={warning}>
                <Warning className={styles.warning}></Warning>
              </Tooltip>}
            </div>
          </div>
          {allActions}
        </div>
      </div>
    </Floatable>
  </Box>
}))

Plot.propTypes = {
  data: PropTypes.any, // Plotly.js data object: Set to null to show Placeholder, set to false to show NoData.
  layout: PropTypes.object, // Plotly.js layout object. Notice that all the fields that are set in this object will be updated.
  config: PropTypes.object, // Plotly.js config object
  capture: PropTypes.object, // Capture settings
  fixedMargins: PropTypes.bool, // When enabled, the margins are not updated upon pan/scroll. This makes the interaction smoother.
  className: PropTypes.string,
  classes: PropTypes.object,
  placeHolderStyle: PropTypes.string, // The CSS class to apply for the Placeholder component.
  noDataStyle: PropTypes.string, // The CSS class to apply for the NoData component.
  floatTitle: PropTypes.string, // The title of the plot shown in floating mode
  metaInfoLink: PropTypes.string, // A link to the metainfo definition
  warning: PropTypes.string, // Optional message displayed under a warning icon on the top-left corner of the element.
  onRelayout: PropTypes.func,
  onRedraw: PropTypes.func,
  onRelayouting: PropTypes.func,
  onHover: PropTypes.func,
  onReset: PropTypes.func,
  onSelected: PropTypes.func,
  onSelecting: PropTypes.func,
  onDeselect: PropTypes.func,
  onClick: PropTypes.func,
  disableDefaultActions: PropTypes.bool, // Whether to show the default plot actions.
  actions: PropTypes.node, // Additional actions.
  onWebglContextLost: PropTypes.func,
  throttleResize: PropTypes.bool, // Whether to throttle resize using requestAnimationFrame
  /**
   * A RxJS Subject for efficient, non-persistent, layout changes that bypass
   * rendering of the component. Should send messages that contain the new
   * layout object.
  */
  layoutSubject: PropTypes.any,
  'data-testid': PropTypes.string,
  showOptions: PropTypes.bool,
  onClickOptions: PropTypes.func
}
Plot.defaultProps = {
  floatTitle: '',
  fixedMargins: true,
  throttleResize: false
}

export default Plot
