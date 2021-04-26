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
import React, {useEffect, useState, useMemo} from 'react'
import { useRecoilValue } from 'recoil'
import PropTypes from 'prop-types'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import clsx from 'clsx'
import {
  Box
} from '@material-ui/core'
import Plot from '../visualization/Plot'
import { convertSI, convertSILabel, mergeObjects } from '../../utils'
import { withErrorHandler } from '../ErrorHandler'
import { normalizationWarning } from '../../config'

const useStyles = makeStyles({
  root: {}
})

function DOS({data, layout, aspectRatio, className, classes, placeholderStyle, unitsState, type, ...other}) {
  // Merge custom layout with default layout
  const units = useRecoilValue(unitsState)
  const initialLayout = useMemo(() => {
    let defaultLayout = {
      yaxis: {
        title: {
          text: `Energy (${convertSILabel('joule', units)})`
        },
        zeroline: false
      },
      xaxis: {
        showexponent: 'first',
        zeroline: false,
        autorange: false
      },
      showlegend: true,
      legend: {
        x: 1,
        xanchor: 'right',
        yanchor: 'bottom',
        y: 0
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units])

  const [finalData, setFinalData] = useState(data)
  const [finalLayout, setFinalLayout] = useState(initialLayout)
  const styles = useStyles(classes)
  const theme = useTheme()

  // Determine if data is normalized
  let normalizedToHOE
  if (data) {
    if (type === 'vibrational') {
      normalizedToHOE = true
    } else {
      if (data.dos_energies_normalized === undefined) {
        normalizedToHOE = false
      } else {
        normalizedToHOE = true
      }
    }
  }

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff as a side effect, the first render containing
  // the placeholders etc. can be done as fast as possible.
  useEffect(() => {
    if (!data) {
      return
    }

    // Determine quantity names
    const norm = type === 'vibrational' ? '' : normalizedToHOE ? '_normalized' : ''
    const energyName = 'dos_energies' + norm
    const valueName = 'dos_values' + norm

    // Convert units and determine range
    let mins = []
    let maxes = []
    let nChannels = data[valueName].length
    const energies = convertSI(data[energyName], 'joule', units, false)
    const values1 = convertSI(data[valueName][0], '1/joule/meter^3', units, false)
    let values2
    mins.push(Math.min(...values1))
    maxes.push(Math.max(...values1))
    if (nChannels === 2) {
      values2 = convertSI(data[valueName][1], '1/joule/meter^3', units, false)
      mins.push(Math.min(...values2))
      maxes.push(Math.max(...values2))
    }
    let range = [Math.min(...mins), Math.max(...maxes)]

    // Create the final data that will be plotted.
    const plotData = []
    if (nChannels === 2) {
      plotData.push(
        {
          x: values2,
          y: energies,
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
    plotData.push(
      {
        x: values1,
        y: energies,
        type: 'scatter',
        mode: 'lines',
        showlegend: false,
        line: {
          color: theme.palette.primary.main,
          width: 2
        }
      }
    )

    // Normalization line
    if (type !== 'vibrational' && normalizedToHOE) {
      plotData.push({
        x: range,
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

    // Compute the final layout
    let computedLayout = mergeObjects(
      {
        xaxis: {
          title: {
            text: convertSILabel(norm ? 'states/joule/meter^3/atom' : 'states/joule/cell', units)
          },
          range: range
        }
      },
      initialLayout
    )
    setFinalData(plotData)
    setFinalLayout(computedLayout)
  }, [data, units, initialLayout, normalizedToHOE, theme, type])

  return (
    <Box className={clsx(styles.root, className)}>
      <Plot
        data={finalData}
        layout={finalLayout}
        aspectRatio={aspectRatio}
        floatTitle="Density of states"
        warning={normalizedToHOE === false ? normalizationWarning : null}
        {...other}
      >
      </Plot>
    </Box>
  )
}

DOS.propTypes = {
  data: PropTypes.any, // section_dos
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  placeholderStyle: PropTypes.string,
  noDataStyle: PropTypes.string,
  unitsState: PropTypes.object, // Recoil atom containing the unit configuration
  type: PropTypes.string // Type of band structure: electronic or vibrational
}
DOS.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler(DOS, 'Could not load density of states.')
