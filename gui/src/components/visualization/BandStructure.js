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
import React, {useState, useEffect, useMemo} from 'react'
import PropTypes from 'prop-types'
import { isFinite } from 'lodash'
import { useTheme } from '@material-ui/core/styles'
import Plot from '../plotting/Plot'
import { add, distance, mergeObjects } from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { withErrorHandler } from '../ErrorHandler'
import { msgNormalizationWarning } from '../../config'
import { getLineStyles } from '../plotting/common'

const BandStructure = React.memo(({
  data,
  layout,
  className,
  units,
  type,
  ...other
}) => {
  const [finalData, setFinalData] = useState(!data ? data : undefined)
  const [pathSegments, setPathSegments] = useState(undefined)
  const energyUnit = useMemo(() => new Unit('joule'), [])
  const theme = useTheme()

  // Determine the energy reference and normalization
  const energyHighestOccupiedNormalized = useMemo(() => {
    if (!data) {
      return [[undefined, undefined]]
    }
    if (type === 'vibrational') {
      return [[0, true]]
    } else {
      return data.map(d => {
        if (!isFinite(d.energy_highest_occupied)) {
          return [0, false]
        } else {
          return [new Quantity(d.energy_highest_occupied, energyUnit).toSystem(units).value(), true]
        }
      })
    }
  }, [data, energyUnit, units, type])

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff as a side effect, the first render containing
  // the placeholders etc. can be done as fast as possible.
  useEffect(() => {
    if (!data) {
      setFinalData(data)
      return
    }

    const energyName = 'energies'
    const kpointName = 'kpoints'
    const plotData = []
    const channels = data.map(d => d.segment[0][energyName].length)
    const lines = getLineStyles(channels.reduce((a, b) => a + b, 0), null, channels.length > 0 ? channels[0] : null).values()
    data.forEach((d, n) => {
      const nChannels = d.segment[0][energyName].length
      const nBands = d.segment[0][energyName][0][0].length

      // Calculate distances in k-space if missing. These distances in k-space
      // define the plot x-axis spacing.
      const tempSegments = []
      if (d.segment[0].k_path_distances === undefined) {
        let length = 0
        for (const segment of d.segment) {
          const k_path_distances = []
          const nKPoints = segment[energyName][0].length
          const start = segment[kpointName][0]
          const end = segment[kpointName].slice(-1)[0]
          const segmentLength = distance(start, end)
          for (let iKPoint = 0; iKPoint < nKPoints; ++iKPoint) {
            const kPoint = segment[kpointName][iKPoint]
            const dist = distance(start, kPoint)
            k_path_distances.push(length + dist)
          }
          length += segmentLength
          segment.k_path_distances = k_path_distances
          tempSegments.push(k_path_distances)
        }
      } else {
        for (const segment of d.segment) {
          tempSegments.push(segment.k_path_distances)
        }
      }
      setPathSegments(tempSegments)

      // Path
      let path = []
      for (const segment of d.segment) {
        path = path.concat(segment.k_path_distances)
        tempSegments.push(segment.k_path_distances)
      }

      //  we assign first for spin up
      const lineUp = lines.next().value
      // Second spin channel
      if (nChannels === 2) {
        const line = lines.next().value
        const bands = []
        for (let iBand = 0; iBand < nBands; ++iBand) {
          bands.push([])
        }
        for (const segment of d.segment) {
          for (let iBand = 0; iBand < nBands; ++iBand) {
            const nKPoints = segment[energyName][1].length
            for (let iKPoint = 0; iKPoint < nKPoints; ++iKPoint) {
              bands[iBand].push(segment[energyName][1][iKPoint][iBand])
            }
          }
        }

        // Create plot data entry for each band
        for (let band of bands) {
          band = new Quantity(band, energyUnit).toSystem(units).value()
          if (energyHighestOccupiedNormalized[n][0] !== 0) {
            band = add(band, -energyHighestOccupiedNormalized[n][0])
          }
          plotData.push(
            {
              x: path,
              y: band,
              type: 'scatter',
              mode: 'lines',
              showlegend: false,
              legendgroup: d.name,
              line: line
            }
          )
        }
      }

      // First spin channel
      const bands = []
      for (let iBand = 0; iBand < nBands; ++iBand) {
        bands.push([])
      }
      for (const segment of d.segment) {
        for (let iBand = 0; iBand < nBands; ++iBand) {
          const nKPoints = segment[energyName][0].length
          for (let iKPoint = 0; iKPoint < nKPoints; ++iKPoint) {
            bands[iBand].push(segment[energyName][0][iKPoint][iBand])
          }
        }
        path = path.concat(segment.k_path_distances)
      }

      // Create plot data entry for each band
      for (let [nBand, band] of bands.entries()) {
        band = new Quantity(band, energyUnit).toSystem(units).value()
        if (energyHighestOccupiedNormalized[n][0] !== 0) {
          band = add(band, -energyHighestOccupiedNormalized[n][0])
        }
        plotData.push(
          {
            x: path,
            y: band,
            name: d.name,
            type: 'scatter',
            mode: 'lines',
            showlegend: nBand === 0 && d.name !== undefined,
            legendgroup: d.name,
            line: lineUp
          }
        )
      }

      // Normalization line
      if (type !== 'vibrational' && energyHighestOccupiedNormalized[0][1] && n === data.length - 1) {
        plotData.push({
          x: [path[0], path[path.length - 1]],
          y: [0, 0],
          name: 'Highest occupied',
          showlegend: true,
          type: 'line',
          mode: 'lines',
          line: {
            color: '#000',
            width: 1
          }
        })
      }
    })

    setFinalData(plotData)
  }, [data, theme, units, energyUnit, type, energyHighestOccupiedNormalized])

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    const defaultLayout = {
      xaxis: {
        tickangle: 0,
        tickfont: {
          size: 14
        },
        zeroline: false
      },
      yaxis: {
        title: {
          text: `Energy (${energyUnit.toSystem(units).label()})`
        },
        zeroline: type === 'vibrational'
      },
      title: {
        text: {
          text: 'Band structure'
        }
      },
      showlegend: true,
      legend: {
        x: 1,
        y: 0,
        xanchor: 'right',
        yanchor: 'bottom',
        font: {
          size: 13
        }
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, type, energyUnit, units])

  // Compute layout that depends on data.
  const computedLayout = useMemo(() => {
    if (data === undefined || data === false || pathSegments === undefined) {
      return {}
    }
    // Set new layout that contains the segment labels
    const labels = []
    const labelKPoints = []
    for (let iSegment = 0; iSegment < data[0].segment.length; ++iSegment) {
      const segment = data[0].segment[iSegment]
      const startLabel = segment.endpoints_labels
        ? segment.endpoints_labels[0]
        : ''
      if (iSegment === 0) {
        labels.push(startLabel)
        labelKPoints.push(pathSegments[iSegment][0])
      } else {
        const prevLabel = labels[labels.length - 1]
        if (prevLabel !== startLabel) {
          labels[labels.length - 1] = `${prevLabel}|${startLabel}`
        }
      }
      const endLabel = segment.endpoints_labels
        ? segment.endpoints_labels[1]
        : ''
      labels.push(endLabel)
      labelKPoints.push(pathSegments[iSegment].slice(-1)[0])
    }

    const shapes = []
    for (let iShape = 1; iShape < labelKPoints.length - 1; ++iShape) {
      const labelKPoint = labelKPoints[iShape]
      shapes.push({
        type: 'line',
        x0: labelKPoint,
        y0: 0,
        x1: labelKPoint,
        y1: 1,
        yref: 'paper',
        line: {
          color: '#999',
          width: 1,
          dash: '10px,10px'
        }
      })
    }
    const ticks = {
      shapes: shapes,
      xaxis: {
        tickmode: 'array',
        tickvals: labelKPoints,
        ticktext: labels
      }
    }
    return ticks
  }, [data, pathSegments])

  // Merge the given layout and layout computed from data
  const finalLayout = useMemo(() => {
    return mergeObjects(computedLayout, tmpLayout)
  }, [computedLayout, tmpLayout])

  return <Plot
    data={finalData}
    layout={finalLayout}
    floatTitle={'Band structure'}
    warning={energyHighestOccupiedNormalized[0][1] ? null : msgNormalizationWarning}
    metaInfoLink={data ? data[0]?.m_path : null}
    className={className}
    {...other}
  >
  </Plot>
})

BandStructure.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.bool, // Set to False to show NoData component
    PropTypes.arrayOf(PropTypes.shape({
      segment: PropTypes.array.isRequired, // Array of segments in SI units
      energy_highest_occupied: PropTypes.number, // Highest occupied energy.
      m_path: PropTypes.string // Path of the section containing the data in the Archive
    }))
  ]),
  layout: PropTypes.object,
  className: PropTypes.string,
  units: PropTypes.object, // Contains the unit configuration
  type: PropTypes.string // Type of band structure: electronic or vibrational
}
BandStructure.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler('Could not load band structure.')(BandStructure)
