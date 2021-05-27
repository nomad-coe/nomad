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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import clsx from 'clsx'
import {
  Box
} from '@material-ui/core'
import Plot from '../visualization/Plot'
import { add, distance, mergeObjects } from '../../utils'
import { convertSI } from '../../units'
import { withErrorHandler } from '../ErrorHandler'
import { normalizationWarning } from '../../config'

const useStyles = makeStyles({
  root: {
  }
})

const BandStructure = React.memo(({
  data,
  layout,
  aspectRatio,
  className,
  classes,
  placeholderStyle,
  units,
  type,
  ...other
}) => {
  const [finalData, setFinalData] = useState(data === false ? data : undefined)
  const [pathSegments, setPathSegments] = useState(undefined)
  const [normalizedToHOE, setNormalizedToHOE] = useState(false)

  // Styles
  const style = useStyles(classes)
  const theme = useTheme()

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff as a side effect, the first render containing
  // the placeholders etc. can be done as fast as possible.
  useEffect(() => {
    if (!data) {
      return
    }

    // Determine the energy reference.
    let energyHighestOccupied
    if (type === 'vibrational') {
      energyHighestOccupied = 0
      setNormalizedToHOE(true)
    } else {
      if (!data.energy_highest_occupied === undefined) {
        energyHighestOccupied = 0
        setNormalizedToHOE(false)
      } else {
        energyHighestOccupied = convertSI(data.energy_highest_occupied, 'joule', units, false)
        setNormalizedToHOE(true)
      }
    }
    const energyName = 'band_energies'
    const kpointName = 'band_k_points'

    let plotData = []
    let nChannels = data.segments[0][energyName].length
    let nBands = data.segments[0][energyName][0][0].length

    // Calculate distances in k-space if missing. These distances in k-space
    // define the plot x-axis spacing.
    let tempSegments = []
    if (data.segments[0].k_path_distances === undefined) {
      let length = 0
      for (let segment of data.segments) {
        const k_path_distances = []
        const nKPoints = segment[energyName][0].length
        let start = segment[kpointName][0]
        let end = segment[kpointName].slice(-1)[0]
        let segmentLength = distance(start, end)
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
      for (let segment of data.segments) {
        tempSegments.push(segment.k_path_distances)
      }
    }
    setPathSegments(tempSegments)

    // Path
    let path = []
    for (let segment of data.segments) {
      path = path.concat(segment.k_path_distances)
      tempSegments.push(segment.k_path_distances)
    }

    // Second spin channel
    if (nChannels === 2) {
      let bands = []
      for (let iBand = 0; iBand < nBands; ++iBand) {
        bands.push([])
      }
      for (let segment of data.segments) {
        for (let iBand = 0; iBand < nBands; ++iBand) {
          let nKPoints = segment[energyName][1].length
          for (let iKPoint = 0; iKPoint < nKPoints; ++iKPoint) {
            bands[iBand].push(segment[energyName][1][iKPoint][iBand])
          }
        }
      }

      // Create plot data entry for each band
      for (let band of bands) {
        band = convertSI(band, 'joule', units, false)
        if (energyHighestOccupied !== 0) {
          band = add(band, -energyHighestOccupied)
        }
        plotData.push(
          {
            x: path,
            y: band,
            type: 'scatter',
            mode: 'lines',
            showlegend: false,
            line: {
              color: theme.palette.secondary.main,
              width: 2
            }
          }
        )
      }
    }

    // First spin channel
    let bands = []
    for (let iBand = 0; iBand < nBands; ++iBand) {
      bands.push([])
    }
    for (let segment of data.segments) {
      for (let iBand = 0; iBand < nBands; ++iBand) {
        let nKPoints = segment[energyName][0].length
        for (let iKPoint = 0; iKPoint < nKPoints; ++iKPoint) {
          bands[iBand].push(segment[energyName][0][iKPoint][iBand])
        }
      }
      path = path.concat(segment.k_path_distances)
    }

    // Create plot data entry for each band
    for (let band of bands) {
      band = convertSI(band, 'joule', units, false)
      if (energyHighestOccupied !== 0) {
        band = add(band, -energyHighestOccupied)
      }
      plotData.push(
        {
          x: path,
          y: band,
          type: 'scatter',
          mode: 'lines',
          showlegend: false,
          line: {
            color: theme.palette.primary.main,
            width: 2
          }
        }
      )
    }

    // Normalization line
    if (type !== 'vibrational' && normalizedToHOE) {
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

    setFinalData(plotData)
  }, [data, theme, units, type, normalizedToHOE])

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    let defaultLayout = {
      xaxis: {
        tickangle: 0,
        tickfont: {
          size: 14
        },
        zeroline: false
      },
      yaxis: {
        title: {
          text: 'Energy (eV)'
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
        yanchor: 'bottom'
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, type])

  // Compute layout that depends on data.
  const computedLayout = useMemo(() => {
    if (data === undefined || data === false || pathSegments === undefined) {
      return {}
    }
    // Set new layout that contains the segment labels
    let labels = []
    let labelKPoints = []
    for (let iSegment = 0; iSegment < data.segments.length; ++iSegment) {
      let segment = data.segments[iSegment]
      const startLabel = segment.band_segm_labels ? segment.band_segm_labels[0] : ''
      if (iSegment === 0) {
        // If label is not defined, use empty string
        labels.push(startLabel)
        labelKPoints.push(pathSegments[iSegment][0])
      } else {
        let prevLabel = labels[labels.length - 1]
        if (prevLabel !== startLabel) {
          labels[labels.length - 1] = `${prevLabel}|${startLabel}`
        }
      }
      const endLabel = segment.band_segm_labels ? segment.band_segm_labels[1] : ''
      labels.push(endLabel)
      labelKPoints.push(pathSegments[iSegment].slice(-1)[0])
    }
    let shapes = []
    for (let iShape = 1; iShape < labelKPoints.length - 1; ++iShape) {
      let labelKPoint = labelKPoints[iShape]
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
    let ticks = {
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

  return (
    <Box className={clsx(style.root, className)}>
      <Plot
        data={finalData}
        layout={finalLayout}
        aspectRatio={aspectRatio}
        floatTitle={'Band structure'}
        warning={normalizedToHOE ? null : normalizationWarning}
        placeholderStyle={placeholderStyle}
        metaInfoLink={data?.m_path}
        {...other}
      >
      </Plot>
    </Box>
  )
})

BandStructure.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.bool, // Set to False to show NoData component
    PropTypes.shape({
      segments: PropTypes.array.isRequired, // Array of section_k_band_segments in SI units
      energy_highest_occupied: PropTypes.number, // Highest occupied energy.
      m_path: PropTypes.string // Path of the section containing the data in the Archive
    })
  ]),
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  placeholderStyle: PropTypes.any,
  units: PropTypes.object, // Contains the unit configuration
  type: PropTypes.string // Type of band structure: electronic or vibrational
}
BandStructure.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler(BandStructure, 'Could not load band structure.')
