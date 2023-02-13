import React, {useMemo} from 'react'
import {Box, useTheme} from '@material-ui/core'
import {Quantity as Q, useUnits} from '../../units'
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
  const units = useUnits()
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
        const quantity = new Q(value, unit).toSystem(units)
        return [quantity.value(), quantity.label()]
      } else {
        return [value, unit]
      }
    }

    /**
     * Recursively splits a path containing slice notation (`start:stop`)
     * @param {String} path The path which should be split at repeating sections
     * @returns {Array} The paths and names of the new paths
     */
    function resolveSlice(path) {
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
              const label = repeat[j][repeatDef?.more?.label_quantity] ?? repeat[j].name
              const [paths, names] = resolveSlice(`${repeatPath}/${j}/${remainPath}`)
              const subSectionName = label ?? `${titleCase(pathArr[i - 1])} ${j}`
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

    const Y = []
    const Names = []
    YPrime.forEach((y) => {
      const [paths, names] = resolveSlice(y)
      Y.push(...paths)
      Names.push(...names)
    })
    const X = []
    XPrime.forEach((x) => {
      X.push(...resolveSlice(x)[0])
    })
    const nLines = Y.length

    const xUnits = []
    const xLabels = []
    const xValuesArray = []
    X.forEach((x) => {
      const [xValues, xUnit] = toUnit(x)
      const xPath = x.split('/')
      const xLabel = titleCase(xPath[xPath.length - 1])
      xUnits.push(xUnit)
      xLabels.push(xLabel)
      xValuesArray.push(xValues)
    })

    const isMultiX = X.length > 1

    const lines = getLineStyles(nLines, theme).map(line => {
      return {type: 'scatter',
        mode: 'lines',
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
      const [yValues, yUnit] = toUnit(y)
      const yPath = y.split('/')
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
