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

function DOS({data, layout, resetLayout, aspectRatio, className, classes, unitsState, ...other}) {
  const [finalData, setFinalData] = useState(undefined)
  const units = useRecoilValue(unitsState)

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    let defaultLayout = {
      yaxis: {
        title: {
          text: `Energy (${convertSILabel('joule', units)})`
        }
      },
      xaxis: {
        showexponent: 'first'
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units])

  // Styles
  const useStyles = makeStyles(
    {
      root: {
      }
    }
  )
  const style = useStyles(classes)
  const theme = useTheme()

  // The plotted data is loaded only after the first render as a side effect to
  // avoid freezing the UI
  useEffect(() => {
    if (data === undefined) {
      return
    }
    const norm = data.dos_energies_normalized === undefined ? '' : '_normalized'
    const energyName = 'dos_energies' + norm
    const valueName = 'dos_values' + norm
    const plotData = []
    if (data !== undefined) {
      let nChannels = data[valueName].length
      let energies = convertSI(data[energyName], 'joule', units, false)
      if (nChannels === 2) {
        plotData.push(
          {
            x: data[valueName][1],
            y: energies,
            type: 'scatter',
            mode: 'lines',
            line: {
              color: theme.palette.secondary.main,
              width: 2
            }
          }
        )
      }
      plotData.push(
        {
          x: data[valueName][0],
          y: energies,
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

  // Compute layout that depends on data.
  const computedLayout = useMemo(() => {
    if (data === undefined) {
      return {}
    }
    const norm = data.dos_energies_normalized !== undefined
    let defaultLayout = {
      xaxis: {
        title: {
          text: norm ? convertSILabel('states/joule/m^3/atom', units) : convertSILabel('states/joule/cell', units)
        }
      }
    }
    return defaultLayout
  }, [data, units])

  // Merge the given layout and layout computed from data
  const finalLayout = useMemo(() => {
    return mergeObjects(computedLayout, tmpLayout)
  }, [computedLayout, tmpLayout])

  return (
    <Box className={clsx(style.root, className)}>
      <Plot
        data={finalData}
        layout={finalLayout}
        resetLayout={resetLayout}
        aspectRatio={aspectRatio}
        floatTitle="Density of states"
        {...other}
      >
      </Plot>
    </Box>
  )
}

DOS.propTypes = {
  data: PropTypes.object, // section_dos
  layout: PropTypes.object,
  resetLayout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  unitsState: PropTypes.object // Recoil atom containing the unit configuration
}

export default withErrorHandler(DOS, 'Could not load density of states.')
