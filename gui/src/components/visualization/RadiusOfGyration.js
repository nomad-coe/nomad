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
import { set, get, isNil, isEqual, isArray, has, capitalize, isPlainObject} from 'lodash'
import { useTheme } from '@material-ui/core/styles'
import Plot from '../plotting/Plot'
import { PropertyItem, PropertySubGrid } from '../entry/properties/PropertyCard'
import { getLocation } from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'
import { ErrorHandler, withErrorHandler } from '../ErrorHandler'
import { getLineStyles } from '../plotting/common'

const timeUnit = new Unit('second')
const distanceUnit = new Unit('meter')
export const rgPath = ['results', 'properties', 'structural', 'radius_of_gyration']
export const rgError = 'Could not load radius of gyration data.'

/**
 * Creates plots for a set of radii of gyration (Rgs). The
 * different types of rgs (molecular, atomic, etc.) get their own subtitle and
 * all Rgs (regardless of label) are plotted in the same graph.
 */
const RadiusOfGyration = React.memo(({
  rg,
  className,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const [finalData, setFinalData] = useState(rg)
  const [finalLayout, setFinalLayout] = useState()
  const {units} = useUnitContext()

  // Check that the data is valid. Otherwise raise an exception.
  assert(isPlainObject(rg), 'Invalid rg data provided.')

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    // Loop through different types
    const tracesAll = {}
    const layoutAll = {}
    for (const [type, values] of Object.entries(rg)) {
      const path = [type]
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
      const traces = [[], []]
      const lineStyles = getLineStyles(values.length, theme)
      const layout = []
      layout.push({
        showlegend: true,
        legend: {
          x: 1,
          y: 1,
          xanchor: 'right'
        },
        xaxis: {
          title: {
            text: `Time (${timeUnit.toSystem(units).label()})`
          },
          autorange: true,
          zeroline: false
        },
        yaxis: {
          title: {
            text: `Rg (${distanceUnit.toSystem(units).label()})`
          },
          autorange: true,
          zeroline: false
        }
      })
      layout.push({
        showlegend: true,
        legend: {
          x: 1,
          y: 1,
          xanchor: 'right'
        },
        yaxis: {
          title: 'Occurrence (%)',
          autorange: true,
          zeroline: false
        },
        xaxis: {
          title: {
            text: `Rg (${distanceUnit.toSystem(units).label()})`
          },
          autorange: true,
          zeroline: false
        }
      })

      for (const value of values) {
        assert(isArray(value?.time), `Radius of gyration was expecting an array of times, but instead got: ${rg?.time}`)
        assert(isArray(value?.value), `Radius of gyration was expecting an array of values, but instead got: ${rg?.value}`)
        const t = new Quantity(value?.time, timeUnit).toSystem(units).value()
        const val = new Quantity(value?.value, distanceUnit).toSystem(units).value()

        traces[0].push({
          x: t,
          y: val,
          name: `${value?.label}`,
          type: 'scatter',
          showlegend: true,
          line: lineStyles[i]
        })

        traces[1].push({
          x: val,
          name: `${value?.label}`,
          type: 'histogram',
          histnorm: 'percent',
          showlegend: true,
          line: lineStyles[i]
        })
        ++i
      }
      set(tracesAll, path, traces)
      set(layoutAll, path, layout)
    }
    setFinalData(tracesAll)
    setFinalLayout(layoutAll)
  }, [rg, units, theme])

  const plots = useMemo(() => {
    const plots = []
    finalData && Object.entries(finalData).forEach(([type, values]) => {
      const title = `${capitalize(type)} radii of gyration`
      plots.push(
        <ErrorHandler message={`Could not load ${title.toLowerCase()}`} key={type}>
          <PropertyItem xs={12} title={title} height="auto" key={type}>
            <PropertySubGrid>
              <PropertyItem xs={12} height="300px">
                <Plot
                  data={values && values?.[0]}
                  layout={finalLayout?.[type]?.[0]}
                  floatTitle={title}
                  fixedMargins={true}
                  className={className}
                  data-testid={`${testID}-${type.toLowerCase()}-0`}
                />
              </PropertyItem>
              <PropertyItem xs={12} height="300px">
                <Plot
                  data={values && values?.[1]}
                  layout={finalLayout?.[type]?.[1]}
                  floatTitle={title}
                  fixedMargins={true}
                  className={className}
                  data-testid={`${testID}-${type.toLowerCase()}-1`}
                />
              </PropertyItem>
            </PropertySubGrid>
          </PropertyItem>
        </ErrorHandler>
      )
    })
    return plots
  }, [className, finalData, finalLayout, testID])

  return plots
})

const rgPlotShape = PropTypes.oneOfType([
  PropTypes.shape({
    time: PropTypes.arrayOf(PropTypes.number),
    value: PropTypes.arrayOf(PropTypes.number)
  }),
  PropTypes.oneOf([false, undefined])
])
const rgLabelShape = PropTypes.arrayOf(rgPlotShape)
const rgShape = PropTypes.objectOf(rgLabelShape)

RadiusOfGyration.propTypes = {
  rg: rgShape,
  methodology: PropTypes.object,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

RadiusOfGyration.defaultProps = {
  'data-testid': 'radius-of-gyration'
}

export default withErrorHandler(rgError)(RadiusOfGyration)

/**
 * Returns the plots for all radii of gyration that are contained in
 * the given archive.
 */
const RadiiOfGyrationRaw = React.memo(({
  index,
  archive
}) => {
  const fromEntry = useCallback((index, archive) => {
    const data = {
      rg: false,
      methodology: false
    }
    const urlPrefix = `${getLocation()}/data`
    const rgIndex = get(index, rgPath)

    // No data
    if (!rgIndex) {
      return data
    }

    // Load the basic information from the index. If clashing methodologies are
    // found, an error is raised for now. In the future several methodologies
    // could be supported by plotting them separately.
    rgIndex.forEach((rg, i) => {
      if (rg.methodology) {
        if (data.methodology && !isEqual(data.methodology, rg.methodology)) {
          throw Error('Multiple methodologies for Rg not yet supported.')
        }
        data.methodology = rg.methodology
      }
      set(data.rg, [rg.kind], undefined)
    })

    // Still loading archive but the amount of plots is already known
    if (!archive) {
      return data
    }

    // Full data ready
    data.rg = {}
    rgIndex.forEach((rg, i) => {
      const kind = rg.kind
      if (!has(data.rg, kind)) {
        data.rg[kind] = []
      }
      const rgArchive = get(archive, rgPath)?.[i]
      data.rg[kind].push({
        ...rg,
        ...rgArchive,
        archiveUrl: `${urlPrefix}/${rgPath.join('/')}:${i}`
      })
    })
    return data
  }, [])
  const data = fromEntry(index, archive)

  return <RadiusOfGyration {...data} />
})

RadiiOfGyrationRaw.propTypes = {
  index: PropTypes.object,
  archive: PropTypes.object
}

export const RadiiOfGyration = withErrorHandler(rgError)(RadiiOfGyrationRaw)
