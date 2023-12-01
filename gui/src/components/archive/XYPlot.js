import React, {useMemo} from 'react'
import {Box, useTheme} from '@material-ui/core'
import {Quantity as Q} from '../units/Quantity'
import { useUnitContext } from '../units/UnitContext'
import {titleCase, resolveInternalRef} from '../../utils'
import {getLineStyles} from '../plotting/common'
import { merge } from 'lodash'
import Plot from '../plotting/Plot'
import PropTypes from 'prop-types'
import {withErrorHandler} from '../ErrorHandler'

class XYPlotError extends Error {
  constructor(message) {
    super(message)
    this.name = 'XYPlotError'
  }
}

const XYPlot = React.memo(function XYPlot({plot, section, sectionDef, title}) {
  const theme = useTheme()
  const {units} = useUnitContext()
  const xAxis = plot.x || plot['x_axis'] || plot['xAxis']
  const yAxis = plot.y || plot['y_axis'] || plot['yAxis']

  const [data, layout] = useMemo(() => {
    if (!sectionDef?._properties) {
      return [undefined, undefined]
    }
    const XPrime = Array.isArray(xAxis) ? xAxis : [xAxis]
    const YPrime = Array.isArray(yAxis) ? yAxis : [yAxis]
    const toUnit = path => {
      const relativePath = '/' + path.replace('./', '')
      const resolvedQuantityDef = resolveInternalRef(relativePath, sectionDef)
      if (resolvedQuantityDef === undefined || resolvedQuantityDef === null) {
        throw new XYPlotError(`Could not resolve the path ${path}`)
      }
      let value
      try {
        value = resolveInternalRef(relativePath, section)
      } catch (err) {
        if (!err.message.startsWith('Path does not exist:')) {
          throw new XYPlotError(err.message)
        }
      }
      if (value === undefined || value === null) {
        // there is not data yet
        return [undefined, undefined]
      }
      const unit = resolvedQuantityDef?.unit
      if (unit) {
        const displayUnit = resolvedQuantityDef?.m_annotations?.eln?.[0]?.defaultDisplayUnit
        const quantity = displayUnit ? new Q(value, unit).to(displayUnit) : new Q(value, unit).toSystem(units)
        return [quantity.value(), quantity.label()]
      } else {
        return [value, unit]
      }
    }

    /**
     * Recursively splits a path containing slice notation (`start:stop`)
     * @param {String} path The path which should be split at repeating sections
     * @param {Boolean} isScalar Whether the quantity at the end of the path is scalar
     * @returns {Array} The paths and names of the new paths
     */
    function resolveSlice(path, isScalar) {
      const pathArr = path.split('/')
      for (let i = 0; i < pathArr.length; i++) {
        const pathSection = pathArr[i]
        if (pathSection.includes(':')) {
          const repeatPath = pathArr.slice(0, i).join('/')
          const repeat = resolveInternalRef('/' + repeatPath.replace('./', ''), section)
          if (repeat === undefined) { return [[], []] }
          // Get start and stop from python slice notation
          let [start, stop, rest] = pathSection.split(':').map((x) => x === '' ? undefined : parseInt(x))
          if (Number.isNaN(start) || Number.isNaN(stop) || rest !== undefined) {
            throw new XYPlotError(`Invalid slice notation: "${pathSection}"`)
          }
          start = start === undefined || start < -repeat.length ? 0 : (start < 0 ? start + repeat.length : start)
          stop = stop === undefined || stop > repeat.length ? repeat.length : (stop < 0 ? stop + repeat.length : stop)
          // If there are initialized sections in this range
          if (stop > start) {
            const remainPath = pathArr.slice(i + 1).join('/')
            const allPaths = []
            const allNames = []
            for (let j = start; j < stop; j++) {
              const repeatDef = resolveInternalRef('/' + repeatPath.replace('./', '') + `/${j}/sub_section`, sectionDef)
              const label = isScalar ? undefined : (repeat[j][repeatDef?.more?.label_quantity] ?? repeat[j].name)
              const [paths, names] = resolveSlice(`${repeatPath}/${j}/${remainPath}`, isScalar)
              const subSectionName = label ?? `${titleCase(pathArr[i - 1])}${isScalar ? '' : ` ${j}`}`
              names.forEach((name) => { allNames.push(`${subSectionName}, ${name}`) })
              allPaths.push(...paths)
            }
            return [allPaths, allNames]
          }
          return [[], []]
        }
      }
      return [[path], [titleCase(pathArr[pathArr.length - 1])]]
    }

    /**
     * Function for getting all the paths for an array of paths that may contain slice
     * notation. The function also returns names to be used as labels if withNames is true.
     * @param {Array} slicedPaths The array of paths that may contatin slice notation
     * @param {Boolean} withNames Whether or not the functiion should return names
     * @returns Array of paths or Array of Array of paths and names
     */
    function getSlicedPaths(slicedPaths, withNames) {
      const allPaths = []
      const Names = []
      slicedPaths.forEach((slicedPath) => {
        const pathRelative = '/' + slicedPath.replace('./', '')
        // Removes the index in the path and resolve the shape of the quantity
        const isScalar = resolveInternalRef(pathRelative.replace(/\/-?\d*:-?\d*/gm, '/0') + '/shape', sectionDef)?.length === 0
        const [paths, names] = resolveSlice(pathRelative, isScalar)
        // If the quanitity is a scalar an Array of paths is added otherwise the paths are added one by one
        allPaths.push(...(isScalar ? [paths] : paths))
        if (withNames) {
          // For scalar quantities the name is taken as the first one
          Names.push(...(isScalar ? [names[0]] : names))
        }
      })
      return withNames ? [allPaths, Names] : allPaths
    }

    const [Y, Names] = getSlicedPaths(YPrime, true)
    const X = getSlicedPaths(XPrime, false)
    const nLines = Y.length

    /**
     * Gets the value array, unit and split path for:
     * 1. A path to an array quanitity
     * or
     * 2. An array of paths to a scalar quantity
     * @param {*} path String path or Array of string paths
     * @returns Array of [value, unit, split path]
     */
    function getValues(path) {
      let Values = []
      let Unit = ''
      let pathArray = ''
      if (Array.isArray(path)) {
        // If path is an Array of paths the quantity is scalar and Values is build from all of them
        path.forEach((Point) => {
          const [Value, PointUnit] = toUnit(Point)
          Values.push(Value)
          Unit = PointUnit
          pathArray = Point.split('/')
        })
      } else {
        [Values, Unit] = toUnit(path)
        pathArray = path.split('/')
      }
      return [Values, Unit, pathArray]
    }

    const xUnits = []
    const xLabels = []
    const xValuesArray = []
    X.forEach((x) => {
      const [xValues, xUnit, xPath] = getValues(x)
      const xLabel = titleCase(xPath[xPath.length - 1])
      xUnits.push(xUnit)
      xLabels.push(xLabel)
      xValuesArray.push(xValues)
    })

    const isMultiX = X.length > 1

    const lines = getLineStyles(nLines, theme).map(line => {
      return {type: 'scatter',
        line: line}
    })
    if (plot.lines) {
      Y.forEach((y, index) => {
        merge(lines[index], plot.lines[index])
      })
    }

    if (isMultiX) {
      if (X.length !== Y.length) {
        throw new XYPlotError('The length of provided x axes and y axes do not match!')
      }
      if (xUnits.some(unit => unit !== xUnits[0])) {
        throw new XYPlotError('Different units are provided for x data. Multi xAxis plot is not supported!')
      }
    }

    const data = []
    const yUnits = []
    const yLabels = []
    Y.forEach((y, index) => {
      const [yValues, yUnit, yPath] = getValues(y)
      // For scalar quantities the default mode is set to 'markers'
      lines[index].mode = lines[index].mode ?? (Array.isArray(y) ? 'markers' : 'lines')
      const yLabel = titleCase(yPath[yPath.length - 1])
      const line = {
        name: Names[index],
        x: isMultiX ? xValuesArray[index] : xValuesArray[0],
        y: yValues,
        ...lines[index]
      }
      data.push(line)
      yUnits.push(yUnit)
      yLabels.push(yLabel)
    })

    const getColor = index => {
      const line = lines[index]
      if ('mode' in line) {
        if (line.mode === 'lines') {
          return {color: line.line?.color}
        } else if (line.mode === 'markers') {
          return {color: line.marker?.color}
        }
      }
      return {color: '#000000'}
    }

    const sameUnit = !yUnits.some(unit => unit !== yUnits[0])

    const layout = {
      xaxis: {
        title: xUnits[0] ? `${xLabels[0]} (${xUnits[0]})` : xLabels[0],
        fixedrange: false
      },
      yaxis: {
        title: sameUnit ? (yUnits[0] ? `${titleCase(title)} (${yUnits[0]})` : titleCase(title)) : (yUnits[0] ? `${yLabels[0]} (${yUnits[0]})` : yLabels[0]),
        titlefont: !sameUnit && nLines > 1 ? getColor(0) : undefined,
        tickfont: !sameUnit && nLines > 1 ? getColor(0) : undefined
      },
      showlegend: sameUnit && nLines > 1,
      legend: {
        x: 1,
        y: 1,
        xanchor: 'right'
      }
    }

    if (!sameUnit) {
      Y.forEach((y, index) => {
        const color = getColor(index)
        if (index > 0) {
          layout[`yaxis${index + 1}`] = {
            title: yUnits[index] ? `${yLabels[index]} (${yUnits[index]})` : yLabels[index],
            anchor: index === 1 ? 'x' : 'free',
            overlaying: 'y',
            side: index % 2 === 0 ? 'left' : 'right',
            titlefont: nLines > 1 ? color : undefined,
            tickfont: nLines > 1 ? color : undefined,
            position: index % 2 === 0 ? 0.1 * index : 1.0 - 0.1 * (index - 1)
          }
          data[index]['yaxis'] = `y${index + 1}`
        }
      })
    }

    if (plot.layout) {
      merge(layout, plot.layout)
    }
    return [data, layout]
  }, [plot.layout, plot.lines, xAxis, yAxis, section, sectionDef, theme, title, units])

  return <Box minWidth={500} height={500}>
    <Plot
      data={data}
      layout={layout}
      floatTitle={title}
      fixedMargins={true}
      config={plot.config}
    />
  </Box>
})
XYPlot.propTypes = {
  plot: PropTypes.object.isRequired,
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  title: PropTypes.string
}

export default withErrorHandler(
  (error) => error.name === 'XYPlotError'
    ? error.message
    : 'Could not load plot due to an unexpected error.'
)(XYPlot)
