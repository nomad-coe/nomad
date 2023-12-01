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
import React, {useState, useEffect, useMemo, useCallback} from 'react'
import PropTypes from 'prop-types'
import assert from 'assert'
import { set, get, isNil, isArray, has, capitalize, isPlainObject} from 'lodash'
import { useTheme } from '@material-ui/core/styles'
import Plot from '../plotting/Plot'
import { PropertyItem, PropertySubGrid } from '../entry/properties/PropertyCard'
import { getLocation } from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'
import { ErrorHandler, withErrorHandler } from '../ErrorHandler'
import { getLineStyles } from '../plotting/common'

const distanceUnit = new Unit('meter')
export const rdfPath = ['results', 'properties', 'structural', 'radial_distribution_function']
export const rdfError = 'Could not load radial distribution function data.'

/**
 * Creates plots for a set of radial distribution functions (RDFs). The
 * different types of rdfs (molecular, atomic, etc.) get their own subtitle and
 * all RDFs with the same label are plotted in the same graph.
 */
const RadialDistributionFunction = React.memo(({
  rdf,
  className,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const [finalData, setFinalData] = useState(rdf)
  const [finalLayout, setFinalLayout] = useState()
  const {units} = useUnitContext()

  // Check that the data is valid. Otherwise raise an exception.
  assert(isPlainObject(rdf), 'Invalid rdf data provided.')

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    // Loop through different types
    const tracesAll = {}
    const layoutAll = {}
    for (const [type, labels] of Object.entries(rdf)) {
      for (const [label, values] of Object.entries(labels)) {
        const path = [type, label]
        if (isNil(values)) {
          set(tracesAll, path, undefined)
          set(layoutAll, path, undefined)
          continue
        }
        if (values === false) {
          set(tracesAll, path, false)
          continue
        }
        let i = 0
        const traces = []
        const lineStyles = getLineStyles(values.length, theme)
        const layout = {
          showlegend: true,
          legend: {
            x: 1,
            y: 1,
            xanchor: 'right'
          },
          xaxis: {
            title: {
              text: `${label} distance (${distanceUnit.toSystem(units).label()})`
            },
            zeroline: false
          },
          yaxis: {
            title: {
              text: `g(r)`
            },
            zeroline: false
          }
        }
        for (const value of values) {
          assert(isArray(value?.bins), `Radial distribution function was expecting an array of bins, but instead got: ${rdf?.bins}`)
          assert(isArray(value?.value), `Radial distribution function was expecting an array of values, but instead got: ${rdf?.value}`)
          const r = new Quantity(value?.bins, distanceUnit).toSystem(units).value()
          traces.push({
            x: r,
            y: value.value,
            name: `Frames: ${value.frame_start}-${value.frame_end}`,
            type: 'scatter',
            showlegend: true,
            line: lineStyles[i]
          })
          ++i
        }
        set(tracesAll, path, traces)
        set(layoutAll, path, layout)
      }
    }
    setFinalData(tracesAll)
    setFinalLayout(layoutAll)
  }, [rdf, units, theme])

  const plots = useMemo(() => {
    const plots = []
    finalData && Object.entries(finalData).forEach(([type, labels]) => {
      const title = `${capitalize(type)} radial distribution functions`
      plots.push(
        <ErrorHandler message={`Could not load ${title.toLowerCase()}`} key={type}>
          <PropertyItem xs={12} title={title} height="auto" key={type}>
            <PropertySubGrid>
              {Object.keys(labels).map((label, i) => <PropertyItem
                xs={12}
                key={`${type}_${label}_${i}`}
                height="300px"
              ><Plot
                data={finalData?.[type]?.[label]}
                layout={finalLayout?.[type]?.[label]}
                floatTitle={title}
                fixedMargins={true}
                className={className}
                data-testid={`${testID}-${type.toLowerCase()}-${label.toLowerCase()}`}
              /></PropertyItem>)}
            </PropertySubGrid>
          </PropertyItem>
        </ErrorHandler>
      )
    })
    return plots
  }, [className, finalData, finalLayout, testID])

  return plots
})

const rdfPlotShape = PropTypes.oneOfType([
  PropTypes.shape({
    bins: PropTypes.arrayOf(PropTypes.number),
    value: PropTypes.arrayOf(PropTypes.number)
  }),
  PropTypes.oneOf([false, undefined])
])
const rdfLabelShape = PropTypes.objectOf(PropTypes.arrayOf(rdfPlotShape))
const rdfShape = PropTypes.objectOf(rdfLabelShape)

RadialDistributionFunction.propTypes = {
  rdf: rdfShape,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

RadialDistributionFunction.defaultProps = {
  'data-testid': 'radial-distribution-function'
}

export default withErrorHandler(rdfError)(RadialDistributionFunction)

/**
 * Returns the plots for all radial distribution functions that are contained in
 * the given archive.
 */
const RadialDistributionFunctionsRaw = React.memo(({
  index,
  archive
}) => {
  const fromEntry = useCallback((index, archive) => {
    let rdfAll = {}
    const urlPrefix = `${getLocation()}/data`
    const rdfIndex = get(index, rdfPath)

    // No data
    if (!rdfIndex) return rdfAll

    // Load the basic information from the index in order to show the plot
    // layout with placeholders
    rdfIndex.forEach((rdf, i) => {
      set(rdfAll, [rdf.type, rdf.label], undefined)
    })

    // Loading archive
    if (!archive) return rdfAll

    // Full data ready
    rdfAll = {}
    rdfIndex.forEach((rdf, i) => {
      const type = rdf.type
      const label = rdf.label
      if (!has(rdfAll, type)) rdfAll[type] = {}
      if (!has(rdfAll[type], label)) rdfAll[type][label] = []
      const rdfArchive = get(archive, rdfPath)?.[i]
      rdfAll[type][label].push({
        ...rdf,
        ...rdfArchive,
        archiveUrl: `${urlPrefix}/${rdfPath.join('/')}:${i}`
      })
    })
    return rdfAll
  }, [])
  const rdf = fromEntry(index, archive)

  return <RadialDistributionFunction rdf={rdf} />
})

RadialDistributionFunctionsRaw.propTypes = {
  index: PropTypes.object,
  archive: PropTypes.object
}

export const RadialDistributionFunctions = withErrorHandler(rdfError)(RadialDistributionFunctionsRaw)
