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
import React, { useState, useLayoutEffect, useMemo, useCallback, useRef } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { cloneDeep, isNil, range } from 'lodash'
import { Responsive as ResponsiveGridLayout } from "react-grid-layout"
import { useResizeDetector } from 'react-resize-detector'
import { Paper } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import { WidgetScatterPlot } from './WidgetScatterPlot'
import { WidgetPeriodicTable } from './WidgetPeriodicTable'
import { WidgetHistogram } from './WidgetHistogram'
import { WidgetTerms } from './WidgetTerms'
import { useSearchContext } from '../SearchContext'

// Breakpoints are designed so that the periodic table is always properly
// visible for each breakpoint/column combination
const periodicTableMinWidth = 710
const periodicTableCols = 12
const cols = {xxl: 36, xl: 30, lg: 24, md: 18, sm: 12}
const breakpoints = Object.fromEntries(Object.entries(cols).map(([size, nCols]) => {
  return [size, periodicTableMinWidth / periodicTableCols * nCols]
}))
const sizes = Object.entries(breakpoints)
const margin = [0, 0]

/**
 * Used for calculating the number of columns based on the current width.
 * @param {number} width The grid width
 * @returns The number of columns.
 */
function getNCols(width) {
  let breakpoint
  for (const [key, size] of sizes) {
    breakpoint = key
    if (width >= size) break
  }
  return cols[breakpoint]
}

/**
 * This class is used to build a layout where items are added one by one and
 * adding a new item does not modify the existing layout. Used to replace the
 * initial layout created by 'react-grid-layout'.
 */
class Layout {
  constructor(items, nCols) {
    this.items = []
    this.nCols = nCols
    for (const item of items) {
      this.add(item)
    }
  }

  /**
   * Used to add an item to the layout. Notice that the order in which the items
   * are added matters.
   * @param {object} item Object containing at least a width and a height
   */
  add(item) {
    // If the item already has a location, skip the automatization.
    if (!isNil(item.x) && !isNil(item.y)) {
      this.items.push(item)
    // If no location specified, place it in the best available location.
    } else {
      // Get all suitable locations. A location is suitable if the item fits in
      // there.
      const locations = this.getAvailableLocations()
      const locationsFit = locations
        .filter((loc) => item.w <= loc.w && item.h <= loc.h)

      // Choose the best one. Currently the topmost location is the best.
      locationsFit.sort((a, b) => a.y - b.y)
      const bestLocation = locationsFit[0]

      // Add item to best location
      this.items.push({...item, x: bestLocation.x, y: bestLocation.y})
    }
  }

  /**
   * Updates the list of available locations. Available locations are
   * represented by a rectangle that has a top-left coordinate, width and
   * height. Notice that this implementation is not very optimal as it tries out
   * every point in the grid.
   */
  getAvailableLocations() {
    const locations = []
    if (this.items.length > 0) {
      const bottom = Math.max(...this.items.map((item) => item.y + item.h))
      const epsilon = 1e-4
      for (const i of range(this.nCols - 1)) {
        for (const j of range(bottom - 1)) {
          const point = {x: i, y: j}
          // See how much horizontal space there is
          let xMax = point.x
          while (!this.withinItem({x: xMax + epsilon, y: point.y + epsilon}) && xMax < this.nCols) {
            xMax += 1
          }
          // See how much vertical space there is
          let yMax = point.y
          while (!this.withinItem({x: point.x + epsilon, y: yMax + epsilon}) && yMax < bottom) {
            yMax += 1
          }
          if (yMax === bottom) yMax = Infinity
          const width = xMax - point.x
          const height = yMax - point.y
          if (width && height) {
            locations.push({...point, w: width, h: height})
          }
        }
      }
      // Add a location at the bottom row that can be used for items that take up
      // all columns.
      locations.push({x: 0, y: bottom, w: this.nCols, h: Infinity})
    } else {
      locations.push({x: 0, y: 0, w: this.nCols, h: Infinity})
    }
    return locations
  }

  /**
   * Checks if the given coordinate is within any of the items. Idea from:
   * https://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
   *
   * @param {object} point Object containin an x and a y coordinate.
   * @returns Whether the coordinate is within an item.
   */
  withinItem(point) {
    function dot(a, b) {
      return a.x * b.x + a.y * b.y
    }
    for (const item of this.items) {
      const AB = {x: item.w, y: 0}
      const AM = {x: point.x - item.x, y: point.y - item.y}
      const BC = {x: 0, y: item.h}
      const BM = {x: point.x - (item.x + item.w), y: point.y - item.y}
      const ABAB = dot(AB, AB)
      const ABAM = dot(AB, AM)
      const BCBC = dot(BC, BC)
      const BCBM = dot(BC, BM)
      const within = (ABAB >= ABAM) && (ABAM >= 0) && (BCBC >= BCBM) && (BCBM >= 0)
      if (within) return true
    }
    return false
  }

  /**
   * Returns the calculated layout.
   * @returns Object containing a list of items.
   */
  getLayout() {
    return this.items
  }
}

/**
 * A component that is used to display all of the widgets in an interactive
 * grid. The items can be moved around and resized.
 */
const useStyles = makeStyles(theme => {
  return {
    root: {
      position: 'relative',
      marginLeft: theme.spacing(-1),
      marginRight: theme.spacing(-1)
    },
    component: {
      position: 'absolute',
      top: theme.spacing(1),
      bottom: theme.spacing(1.5),
      left: theme.spacing(1.5),
      right: theme.spacing(1.5),
      height: 'unset',
      width: 'unset'
    },
    containerOuter: {
      position: 'relative'
    },
    containerInner: {
      position: 'absolute',
      top: theme.spacing(1),
      bottom: theme.spacing(1),
      left: theme.spacing(1),
      right: theme.spacing(1),
      height: 'unset',
      width: 'unset'
    }
  }
})
const WidgetGrid = React.memo(({
  className,
  classes
}) => {
  const styles = useStyles(classes)
  const { useWidgetsState } = useSearchContext()
  const { ref, width } = useResizeDetector()
  const [widgets, setWidgets] = useWidgetsState()
  const [validWidth, setValidWidth] = useState()
  const nCols = getNCols(validWidth)
  const firstLayout = useRef(true)

  // Create the layouts
  const layout = useMemo(() => {
    if (!nCols) return {}

    const layouts = {}
    for (const breakpoint of Object.keys(breakpoints)) {
      // This is the layout in the format as react-grid-layout would read it. x:
      // Infinity means that we want to place the item at the very end.
      // Add widgets
      let i = 0
      let layout = Object.entries(widgets)
        .filter(([id, value]) => value?.visible)
        .map(([id, value]) => {
          const layout = value.layout?.[breakpoint]
          const config = {
            i: id,
            x: isNil(layout?.x) ? Infinity : layout?.x,
            y: isNil(layout?.y) ? 0 : layout?.y,
            w: isNil(layout?.w) ? 3 : layout?.w,
            h: isNil(layout?.h) ? 3 : layout?.h,
            minW: isNil(layout?.minW) ? 3 : layout?.minW,
            minH: isNil(layout?.minH) ? 3 : layout?.minH,
            index: i
          }
          ++i
          return config
        })

      // Calculate a sane initial layout using the custom Layout class. The
      // initial layout as calculated by react-grid-layout is not always very
      // great.
      if (firstLayout.current) {
        layout = new Layout(layout.sort((a, b) => a.index - b.index), nCols).getLayout()
        firstLayout.current = false
      }
      layouts[breakpoint] = layout
    }
    return layouts
  }, [widgets, nCols])

  // The grid children are memoized as instructed in the docs of
  // 'react-grid-layout' for performance.
  const children = useMemo(() => {
    return Object.entries(widgets)
        .filter(([id, value]) => value?.visible)
        .sort((a, b) => b[1].index - a[1].index)
        .map(([id, value]) => {
          const Comp = {
            scatterplot: WidgetScatterPlot,
            periodictable: WidgetPeriodicTable,
            histogram: WidgetHistogram,
            terms: WidgetTerms
          }[value.type]
          return <div key={id} className={styles.containerOuter}>
            <Paper className={styles.containerInner}>
              <Comp
                id={id}
                {...value}
                className={styles.component}
              ></Comp>
            </Paper>
          </div>
        })
  }, [widgets, styles])

  // Mount the grid only after the width has been succesfully calculated. Also
  // the latest valid width is stored and used for the actual grid, since upon
  // the changing route, the width of the component goes to zero causing
  // unwanted flickering.
  useLayoutEffect(() => {
    if (width) {
      setValidWidth(width)
    }
  }, [width])

  // The layouts are stored when they are changed. This allows retaining the
  // layout for each breakpoint individually. Notice that we need to feed in the
  // changed layout back to the component so that it works in a controlled
  // manner.
  const handleLayoutChange = useCallback((layout, allLayouts) => {
    setWidgets(old => {
      const newLayout = cloneDeep(old)
      for (const [breakpoint, value] of Object.entries(allLayouts)) {
        for (const layout of value) {
          newLayout[layout.i].layout[breakpoint] = cloneDeep(layout)
        }
      }
      return newLayout
    })
  }, [setWidgets])

  return <div ref={ref} className={clsx(className, styles.root)}>
    { validWidth
      ? <ResponsiveGridLayout
        className="layout"
        layouts={layout}
        compactType="horizontal"
        draggableHandle=".dragHandle"
        autoSize={true}
        onLayoutChange={handleLayoutChange}
        isDraggable={true}
        isResizable={true}
        isBounded={true}
        rowHeight={validWidth / nCols}
        margin={margin}
        breakpoints={breakpoints}
        cols={cols}
        resizeHandle={<ResizeHandle/>}
        width={validWidth}
      >
        {children}
      </ResponsiveGridLayout>
    : null }
  </div>
})

WidgetGrid.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object
}

export default WidgetGrid

/**
 * A custom resize handle component for the grid items.
 */
const useHandleStyles = makeStyles(theme => {
  return {
    root: {
      margin: theme.spacing(1.25)
    }
  }
})
const ResizeHandle = React.forwardRef((props, ref) => {
  const {handleAxis, ...restProps} = props
  const styles = useHandleStyles()
  return <div
    ref={ref}
    className={clsx('react-resizable-handle', `react-resizable-handle-${handleAxis}`, styles.root)}
    style={{width: '30px', height: '30px'}}
    {...restProps}
    >
    </div>
})
ResizeHandle.propTypes = {
  handleAxis: PropTypes.string
}
