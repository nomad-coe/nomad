import React, {useMemo} from 'react'
import {Box} from '@material-ui/core'
import {Quantity as Q} from '../units/Quantity'
import {useUnitContext} from '../units/UnitContext'
import {titleCase, resolveInternalRef} from '../../utils'
import {cloneDeep, merge} from 'lodash'
import Plot from '../plotting/Plot'
import PropTypes from 'prop-types'
import {withErrorHandler} from '../ErrorHandler'

class XYPlotError extends Error {
  constructor(message) {
    super(message)
    this.name = 'XYPlotError'
  }
}

const traverse = (value, callback, parent = null, key = null) => {
  if (value && typeof value === 'object') {
    Object.entries(value).forEach(([k, v]) => {
      traverse(v, callback, value, k)
    })
  } else if (value && Array.isArray(value)) {
    value.forEach((item, index) => {
      traverse(item, callback, value, index)
    })
  } else {
    callback(parent, key, value)
  }
}

const PlotlyFigure = React.memo(function PlotlyFigure({plot, section, sectionDef, title, metaInfoLink}) {
  const {units} = useUnitContext()

  const plotlyGraphObj = useMemo(() => {
    if (!sectionDef?._properties) {
      return [undefined, undefined]
    }

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
    function getSlicedPaths(path, withNames) {
      const slicedPaths = Array.isArray(path) ? path : [path]
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
    const yUnits = []
    const yLabels = []
    const zUnits = []
    const zLabels = []

    const resolveReferences = (data, key, units, labels) => {
      const value = data?.[key]
      if (value && (typeof value === 'string' || value instanceof String)) {
        if (value.startsWith("#")) {
          const [paths, names] = getSlicedPaths(value.slice(1), true)
          const [values, unit, path] = getValues(paths[0])
          const label = names[0] || titleCase(path[path.length - 1])
          data[key] = values
          units.push(unit)
          labels.push(label)
        }
      }
    }

    const plotlyGraphObj = cloneDeep(plot)
    const isDataArray = plotlyGraphObj?.data && Array.isArray(plotlyGraphObj.data)
    const dataArray = isDataArray ? plotlyGraphObj.data : [plotlyGraphObj.data]
    dataArray.forEach((data, index) => {
      resolveReferences(data, 'x', xUnits, xLabels)
      resolveReferences(data, 'y', yUnits, yLabels)
      resolveReferences(data, 'z', zUnits, zLabels)
    })

    traverse(plotlyGraphObj.data, (data, key, value) => {
      if (key === 'color') {
        if (typeof value === 'string' || value instanceof String) {
          if (value.startsWith("#")) {
            let paths, values
            try {
              [paths] = getSlicedPaths(value.slice(1), true);
              [values] = getValues(paths[0])
            } catch (error) {
              if (value && /^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/.test(value)) {
                values = value
              } else {
                throw new XYPlotError(error)
              }
            }
            data[key] = values
          }
        }
      }
    })

    const layout = {}
    if (xLabels.length > 0) {
      layout.xaxis = {
        title: xUnits[0] ? `${xLabels[0]} (${xUnits[0]})` : xLabels[0]
      }
    }
    yLabels.forEach((yLabel, index) => {
      layout[index === 0 ? 'yaxis' : `yaxis${index + 1}`] = {
        title: yUnits[index] ? `${yLabel} (${yUnits[index]})` : yLabel
      }
    })
    zLabels.forEach((zLabel, index) => {
      layout[index === 0 ? 'zaxis' : `zaxis${index + 1}`] = {
        title: zUnits[index] ? `${zLabel} (${zUnits[index]})` : zLabel
      }
    })

    plotlyGraphObj.layout = merge({}, layout, plotlyGraphObj.layout)
    return plotlyGraphObj
  }, [plot, section, sectionDef, units])

  return <Box minWidth={500} height={plotlyGraphObj?.layout?.height || plotlyGraphObj?.layout?.template?.layout?.height || 500}>
    <Plot
      data={plotlyGraphObj.data}
      layout={plotlyGraphObj.layout}
      floatTitle={title}
      fixedMargins={true}
      config={plotlyGraphObj.config}
      metaInfoLink={metaInfoLink}
    />
  </Box>
})
PlotlyFigure.propTypes = {
  plot: PropTypes.object.isRequired,
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  title: PropTypes.string,
  metaInfoLink: PropTypes.string
}

export default withErrorHandler(
  (error) => error.name === 'XYPlotError'
    ? error.message
    : 'Could not load plot due to an unexpected error.'
)(PlotlyFigure)
