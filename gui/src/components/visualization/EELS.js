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
import { useTheme } from '@material-ui/core/styles'
import Plot from '../visualization/Plot'
import { getLineStyles, mergeObjects } from '../../utils'
import { Quantity, Unit } from '../../units'
import { withErrorHandler } from '../ErrorHandler'

const energyUnit = new Unit('joule')

/**
 * Graph for EELS (electron energy loss specroscopy) data.
 */
function EELS({data, layout, className, units, ...other}) {
  const [finalData, setFinalData] = useState(undefined)
  const theme = useTheme()

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    if (data === undefined) {
      return
    }
    const nTraces = data.length
    const defaultLayout = {
      showlegend: nTraces > 1,
      legend: {
        x: 1,
        y: 1,
        xanchor: 'right',
        yanchor: 'top'
      },
      yaxis: {
        title: {
          text: 'Electron count'
        }
      },
      xaxis: {
        showexponent: 'first',
        title: {
          text: `Electron energy loss (${energyUnit.toSystem(units).label()})`
        }
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units, data])

  // The plotted data is loaded only after the first render as a side effect to
  // avoid freezing the UI
  useEffect(() => {
    if (data === undefined) {
      return
    }
    const plotData = []
    const lineStyles = getLineStyles(data.length, theme)
    for (let i = 0; i < data.length; ++i) {
      const trace = data[i]
      const energies = new Quantity(trace.energy, energyUnit).toSystem(units).value()
      plotData.push(
        {
          x: energies,
          y: trace.count,
          name: i,
          type: 'scatter',
          mode: 'lines',
          line: lineStyles[i]
        }
      )
    }
    setFinalData(plotData)
  }, [data, theme, units])

  return <Plot
    data={finalData}
    layout={tmpLayout}
    floatTitle="EELS"
    fixedMargins={true}
    className={className}
    {...other}
  />
}

EELS.propTypes = {
  data: PropTypes.arrayOf(PropTypes.shape({
    count: PropTypes.arrayOf(PropTypes.number),
    energy: PropTypes.arrayOf(PropTypes.number)
  })),
  layout: PropTypes.object,
  className: PropTypes.string,
  units: PropTypes.object
}

export default withErrorHandler(EELS, 'Could not load EELS data.')
