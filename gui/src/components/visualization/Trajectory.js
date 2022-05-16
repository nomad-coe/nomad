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
import React, {useMemo} from 'react'
import PropTypes from 'prop-types'
import { useTheme } from '@material-ui/core/styles'
import { Box } from '@material-ui/core'
import Plot from './Plot'
import { mergeObjects } from '../../utils'
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'
import { SectionTableAutomatic } from '../Quantity'
import { Quantity, Unit } from '../../units'
import { ErrorHandler, withErrorHandler } from '../ErrorHandler'
import { getPlotLayoutVertical, getPlotTracesVertical } from '../plotting/common'

const timeUnit = new Unit('second')
const energyUnit = new Unit('joule')
const temperatureUnit = new Unit('kelvin')
const pressureUnit = new Unit('pascal')

/**
 * Graph for the thermodynamic properties reported during molecular dynamics.
 */
function Trajectory({temperature, pressure, energyPotential, methodology, layout, aspectRatio, className, units, archiveURL, ...other}) {
  const theme = useTheme()

  // Create a configuration for the plot data
  const plots = useMemo(() => [
    {
      data: temperature && {
        x: new Quantity(temperature.time, timeUnit).toSystem(units).value(),
        y: new Quantity(temperature.value, temperatureUnit).toSystem(units).value()
      },
      ytitle: {text: `Temperature (${temperatureUnit.toSystem(units).label()})`}
    },
    {
      data: pressure && {
        x: new Quantity(pressure.time, timeUnit).toSystem(units).value(),
        y: new Quantity(pressure.value, pressureUnit).toSystem(units).value()
      },
      ytitle: {text: `Pressure (${pressureUnit.toSystem(units).label()})`}
    },
    {
      data: energyPotential && {
        x: new Quantity(energyPotential.time, timeUnit).toSystem(units).value(),
        y: new Quantity(energyPotential.value, energyUnit).toSystem(units).value()
      },
      ytitle: {text: `Potential energy (${energyUnit.toSystem(units).label()})`}
    }
  ], [energyPotential, pressure, temperature, units])

  // Create final layout. Notice that layout is created before the actual data
  // is loaded.
  const mergedLayout = useMemo(() => {
    const defaultLayout = getPlotLayoutVertical(
      plots,
      {text: `Time (${timeUnit.toSystem(units).label()})`},
      theme
    )
    return mergeObjects(layout, defaultLayout)
  }, [plots, units, theme, layout])

  const nPlots = mergedLayout.grid.rows

  return <PropertyGrid>
    {nPlots !== 0 &&
      <PropertyItem xs={12} title="Trajectory" height="auto">
        <div style={{height: `${nPlots * 200}px`}}>
          <ErrorHandler message='Could not load trajectory data.'>
            <Plot
              data={getPlotTracesVertical(plots, theme)}
              layout={mergedLayout}
              floatTitle="Trajectory"
              fixedMargins={true}
              className={className}
              data-testid='trajectory'
              metaInfoLink={archiveURL}
              {...other}
            />
          </ErrorHandler>
        </div>
        <Box marginTop={1}>
          <SectionTableAutomatic
            data={methodology}
            prefix="results.properties.thermodynamic.trajectory.methodology"
            columns={4}
          />
        </Box>
      </PropertyItem>
    }
  </PropertyGrid>
}

const dynamicShape = PropTypes.oneOfType([
  PropTypes.bool,
  PropTypes.shape({
    value: PropTypes.arrayOf(PropTypes.number),
    time: PropTypes.arrayOf(PropTypes.number)
  })
])

Trajectory.propTypes = {
  temperature: dynamicShape,
  pressure: dynamicShape,
  energyPotential: dynamicShape,
  methodology: PropTypes.object,
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  units: PropTypes.object,
  archiveURL: PropTypes.string // Path for the data in the archive browser
}

export default withErrorHandler(Trajectory, 'Could not load trajectory data.')
