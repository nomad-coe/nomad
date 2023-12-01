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
import PropTypes from 'prop-types'
import { useTheme } from '@material-ui/core/styles'
import Plot from '../plotting/Plot'
import { mergeObjects } from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { withErrorHandler } from '../ErrorHandler'

const HeatCapacity = React.memo(({
  data,
  layout,
  className,
  units,
  'data-testid': testID,
  ...other
}) => {
  const tempUnit = useMemo(() => new Unit('kelvin'), [])
  const capacityUnit = useMemo(() => new Unit('joule / kelvin'), [])
  const [finalData, setFinalData] = useState(!data ? data : undefined)
  const theme = useTheme()

  // Merge custom layout with default layout
  const finalLayout = useMemo(() => {
    const defaultLayout = {
      xaxis: {
        title: {
          text: `Temperature (${tempUnit.toSystem(units).label()})`
        },
        zeroline: false
      },
      yaxis: {
        title: {
          text: `Heat capacity (${capacityUnit.toSystem(units).label()})`
        },
        zeroline: false
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units, tempUnit, capacityUnit])

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    if (!data) {
      setFinalData(data)
      return
    }

    // Convert units and determine range
    const temperatures = new Quantity(data.temperatures, tempUnit).toSystem(units).value()
    const heatCapacities = new Quantity(data.heat_capacities, capacityUnit).toSystem(units).value()

    // Create the final data that will be plotted.
    const plotData = [{
      x: temperatures,
      y: heatCapacities,
      type: 'scatter',
      mode: 'lines',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]

    setFinalData(plotData)
  }, [data, units, tempUnit, capacityUnit, theme])

  return <Plot
    data={finalData}
    layout={finalLayout}
    floatTitle="Heat capacity"
    metaInfoLink={data?.m_path}
    data-testid={testID}
    className={className}
    {...other}
  />
})

HeatCapacity.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.bool, // Set to false to show NoData component
    PropTypes.shape({
      heat_capacities: PropTypes.array.isRequired,
      temperatures: PropTypes.array.isRequired,
      m_path: PropTypes.string // Path of the section containing the data in the Archive
    })
  ]),
  layout: PropTypes.object,
  className: PropTypes.string,
  units: PropTypes.object, // Contains the unit configuration
  'data-testid': PropTypes.string
}
HeatCapacity.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler('Could not load heat capacity.')(HeatCapacity)
