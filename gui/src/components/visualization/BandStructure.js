import React, {useState, useEffect, useMemo} from 'react'
import { useRecoilValue } from 'recoil'
import PropTypes from 'prop-types'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import clsx from 'clsx'
import {
  Box
} from '@material-ui/core'
import Plot from '../visualization/Plot'
import { convertSI, distance, mergeObjects } from '../../utils'

export default function BandStructure({data, layout, aspectRatio, className, classes, onRelayout, onAfterPlot, onRedraw, onRelayouting, unitsState}) {
  const [finalData, setFinalData] = useState(undefined)
  const [pathSegments, setPathSegments] = useState(undefined)
  const units = useRecoilValue(unitsState)

  // Styles
  const useStyles = makeStyles(
    {
      root: {
      }
    }
  )
  const style = useStyles(classes)
  const theme = useTheme()

  // Determine the final plotted data based on the received data. Will work with
  // normalized and unnormalized data.
  useEffect(() => {
    if (data === undefined) {
      return
    }

    // Determine if data is normalized
    const norm = data.section_k_band_segment_normalized === undefined ? '' : '_normalized'
    const segmentName = 'section_k_band_segment' + norm
    const energyName = 'band_energies' + norm
    const kpointName = 'band_k_points' + norm

    let plotData = []
    let nChannels = data[segmentName][0][energyName].length
    let nBands = data[segmentName][0][energyName][0][0].length

    // Calculate distances if missing
    let tempSegments = []
    if (data[segmentName][0].k_path_distances === undefined) {
      let length = 0
      for (let segment of data[segmentName]) {
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
      for (let segment of data[segmentName]) {
        tempSegments.push(segment.k_path_distances)
      }
    }
    setPathSegments(tempSegments)

    // Path
    let path = []
    for (let segment of data[segmentName]) {
      path = path.concat(segment.k_path_distances)
      tempSegments.push(segment.k_path_distances)
    }

    // Second spin channel
    if (nChannels === 2) {
      let bands = []
      for (let iBand = 0; iBand < nBands; ++iBand) {
        bands.push([])
      }
      for (let segment of data[segmentName]) {
        for (let iBand = 0; iBand < nBands; ++iBand) {
          let nKPoints = segment[energyName][0].length
          for (let iKPoint = 0; iKPoint < nKPoints; ++iKPoint) {
            bands[iBand].push(segment[energyName][0][iKPoint][iBand])
          }
        }
      }

      // Create plot data entry for each band
      for (let band of bands) {
        band = convertSI(band, 'joule', units, false)
        plotData.push(
          {
            x: path,
            y: band,
            type: 'scatter',
            mode: 'lines',
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
    for (let segment of data[segmentName]) {
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
      plotData.push(
        {
          x: path,
          y: band,
          type: 'scatter',
          mode: 'lines',
          line: {
            color: theme.palette.primary.main,
            width: 2
          }
        }
      )
    }

    setFinalData(plotData)
  }, [data, theme.palette.primary.main, theme.palette.secondary.main, units])

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    let defaultLayout = {
      xaxis: {
        tickangle: 0,
        tickfont: {
          size: 14
        }
      },
      yaxis: {
        title: {
          text: 'Energy (eV)'
        }
      }
    }
    mergeObjects(layout, defaultLayout)
    return defaultLayout
  }, [layout])

  // Compute layout that depends on data.
  const computedLayout = useMemo(() => {
    if (data === undefined || pathSegments === undefined) {
      return {}
    }
    // Set new layout that contains the segment labels
    const norm = data.section_k_band_segment_normalized === undefined ? '' : '_normalized'
    const segmentName = 'section_k_band_segment' + norm
    const labelName = 'band_segm_labels' + norm
    let labels = []
    let labelKPoints = []
    for (let iSegment = 0; iSegment < data[segmentName].length; ++iSegment) {
      let segment = data[segmentName][iSegment]
      if (iSegment === 0) {
        // If label is not defined, use empty string
        const startLabel = segment[labelName] ? segment[labelName][0] : ''
        labels.push(startLabel)
        labelKPoints.push(pathSegments[iSegment][0])
      }
      const endLabel = segment[labelName] ? segment[labelName][1] : ''
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
    return mergeObjects(computedLayout, tmpLayout, 'shallow')
  }, [computedLayout, tmpLayout])

  return (
    <Box className={clsx(style.root, className)}>
      <Plot
        data={finalData}
        layout={finalLayout}
        aspectRatio={aspectRatio}
        floatTitle={'Band structure'}
        onRelayout={onRelayout}
        onAfterPlot={onAfterPlot}
        onRedraw={onRedraw}
        onRelayouting={onRelayouting}
      >
      </Plot>
    </Box>
  )
}

BandStructure.propTypes = {
  data: PropTypes.object, // section_band_structure or section_band_structure_normalized
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  onAfterPlot: PropTypes.func,
  onRedraw: PropTypes.func,
  onRelayout: PropTypes.func,
  onRelayouting: PropTypes.func,
  unitsState: PropTypes.object // Recoil atom containing the unit configuration
}
