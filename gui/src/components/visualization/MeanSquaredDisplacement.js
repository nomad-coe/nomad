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
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'
import { getLocation, formatNumber, DType } from '../../utils'
import { Quantity, Unit, useUnits } from '../../units'
import { ErrorHandler, withErrorHandler } from '../ErrorHandler'
import { getLineStyles } from '../plotting/common'

const timeUnit = new Unit('second')
const distanceUnitSquared = new Unit('meter**2')
const diffusionUnit = new Unit('meter**2 / second')
export const msdPath = ['results', 'properties', 'dynamical', 'mean_squared_displacement']
export const msdError = 'Could not load mean squared displacement data.'

/**
 * Creates plots for a set of mean squared displacements (msds). The
 * different types of msds (molecular, atomic, etc.) get their own subtitle and
 * all MSDs (regardless of label) are plotted in the same graph.
 */
const MeanSquaredDisplacement = React.memo(({
  msd,
  className,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const [finalData, setFinalData] = useState(msd)
  const [finalLayout, setFinalLayout] = useState()
  const units = useUnits()

  // Check that the data is valid. Otherwise raise an exception.
  assert(isPlainObject(msd), 'Invalid msd data provided.')

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    // Loop through different types
    const tracesAll = {}
    const layoutAll = {}
    for (const [type, values] of Object.entries(msd)) {
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
            text: `time (${timeUnit.toSystem(units).label()})`
          },
          type: 'log',
          autorange: true,
          zeroline: false
        },
        yaxis: {
          title: {
            text: `msd (${distanceUnitSquared.toSystem(units).label()})`
          },
          type: 'log',
          autorange: true,
          zeroline: false
        }
      }
      for (const value of values) {
        assert(isArray(value?.times), `Mean squared displacent was expecting an array of times, but instead got: ${msd?.times}`)
        assert(isArray(value?.value), `Mean squared displacement was expecting an array of values, but instead got: ${msd?.value}`)
        const t = new Quantity(value?.times, timeUnit).toSystem(units).value()
        const val = new Quantity(value?.value, distanceUnitSquared).toSystem(units).value()
        const diffusionConstant = new Quantity(value?.diffusion_constant_value, diffusionUnit).toSystem(units)
        const D = formatNumber(diffusionConstant.value(), DType.Float, 'scientific', 2)
        const D_label = diffusionConstant.label()
        const diffusionConstantError = new Quantity(value?.diffusion_constant_errors, diffusionUnit)
        const R = formatNumber(diffusionConstantError.value(), DType.Float, 'standard', 2)

        traces.push({
          x: t,
          y: val,
          name: `${value?.label}: D=${D} ${D_label}; R=${R}`,
          type: 'scatter',
          showlegend: true,
          line: lineStyles[i]
        })
        ++i
      }
      set(tracesAll, path, traces)
      set(layoutAll, path, layout)
      // }
    }
    setFinalData(tracesAll)
    setFinalLayout(layoutAll)
  }, [msd, units, theme])

  const plots = useMemo(() => {
    const plots = []
    finalData && Object.entries(finalData).forEach(([type, values]) => {
      const title = `${capitalize(type)} mean squared displacements`
      plots.push(
        <ErrorHandler message={`Could not load ${title.toLowerCase()}`} key={type}>
          <PropertyItem xs={12} title={title} height="500px" key={type}>
            <Plot
              data={finalData?.[type]}
              layout={finalLayout?.[type]}
              floatTitle={title}
              fixedMargins={true}
              className={className}
              data-testid={`${testID}-${type.toLowerCase()}`}
            /></PropertyItem>
        </ErrorHandler>
      )
    })
    return plots
  }, [className, finalData, finalLayout, testID])

  return <PropertyGrid>
    {plots}
  </PropertyGrid>
})

const msdPlotShape = PropTypes.oneOfType([
  PropTypes.oneOf([false, undefined]),
  PropTypes.shape({
    times: PropTypes.arrayOf(PropTypes.number),
    value: PropTypes.arrayOf(PropTypes.number)
  })
])
const msdLabelShape = PropTypes.arrayOf(msdPlotShape)
const msdShape = PropTypes.objectOf(msdLabelShape)

MeanSquaredDisplacement.propTypes = {
  msd: msdShape,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}

MeanSquaredDisplacement.defaultProps = {
  'data-testid': 'mean-squared-displacement'
}

export default withErrorHandler(msdError)(MeanSquaredDisplacement)

/**
 * Returns the plots for all mean squared displacements that are contained in
 * the given archive.
 */
const MeanSquaredDisplacementsRaw = React.memo(({
  index,
  archive
}) => {
  const fromEntry = useCallback((index, archive) => {
    let msdAll = {}
    const urlPrefix = `${getLocation()}/data`
    const msdIndex = get(index, msdPath)

    // No data
    if (!msdIndex) return msdAll

    // Load the basic information from the index in order to show the plot
    // layout with placeholders
    msdIndex.forEach((msd) => {
      set(msdAll, [msd.type], undefined)
    })

    // Loading archive
    if (!archive) return msdAll

    // Full data ready
    msdAll = {}
    msdIndex.forEach((msd, i) => {
      const type = msd.type
      if (!has(msdAll, type)) msdAll[type] = []
      const msdArchive = get(archive, msdPath)?.[i]
      msdAll[type].push({
        ...msd,
        ...msdArchive,
        archiveUrl: `${urlPrefix}/${msdPath.join('/')}:${i}`
      })
    })
    return msdAll
  }, [])
  const msd = fromEntry(index, archive)

  return <MeanSquaredDisplacement msd={msd} />
})

MeanSquaredDisplacementsRaw.propTypes = {
  index: PropTypes.object,
  archive: PropTypes.object
}

export const MeanSquaredDisplacements = withErrorHandler(msdError)(MeanSquaredDisplacementsRaw)
