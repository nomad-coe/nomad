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
import React, { useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { isNil, get } from 'lodash'
import { useTheme, makeStyles } from '@material-ui/core/styles'
import { Box } from '@material-ui/core'
import Plot from '../plotting/Plot'
import { mergeObjects, getLocation } from '../../utils'
import {
  PropertyMethodologyItem,
  PropertyMethodologyList
} from '../entry/properties/PropertyCard'
import { Quantity, Unit, useUnits } from '../../units'
import { withErrorHandler } from '../ErrorHandler'
import { getPlotLayoutVertical, getPlotTracesVertical } from '../plotting/common'

const timeUnit = new Unit('second')
const energyUnit = new Unit('joule')
const temperatureUnit = new Unit('kelvin')
const pressureUnit = new Unit('pascal')
export const trajectoryPath = ['results', 'properties', 'thermodynamic', 'trajectory']
export const trajectoryError = 'Could not load trajectory data.'

/**
 * Plot for the thermodynamic properties reported during molecular dynamics.
 */
const useStyles = makeStyles((theme) => ({
  trajectory: {}
}))
const Trajectory = React.memo(({
  temperature,
  pressure,
  energyPotential,
  methodology,
  layout,
  className,
  classes,
  archiveURL,
  ...other
}) => {
  let nPlots = 0
  if (temperature !== false) ++nPlots
  if (pressure !== false) ++nPlots
  if (energyPotential !== false) ++nPlots
  const styles = useStyles({classes: classes})
  const theme = useTheme()
  const units = useUnits()
  const [finalData, setFinalData] = useState(nPlots === 0 ? false : undefined)
  const [finalLayout, setFinalLayout] = useState()

  // Create a configuration for the plot data
  useEffect(() => {
    // Since the properties are plotted together, if any of them is laoding, show the
    // placeholder
    if ([temperature, pressure, energyPotential].some(isNil)) {
      setFinalData(undefined)
      return
    }
    // If there is no data, show the NoData component.
    if (nPlots === 0) {
      setFinalData(false)
      return
    }
    const traces = [
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
    ]
    const defaultLayout = getPlotLayoutVertical(
      traces,
      {text: `Time (${timeUnit.toSystem(units).label()})`},
      theme
    )
    const mergedLayout = mergeObjects(layout, defaultLayout)
    setFinalData(getPlotTracesVertical(traces, theme))
    setFinalLayout(mergedLayout)
  }, [temperature, pressure, energyPotential, layout, units, theme, nPlots])

  return <Box display="flex" flexDirection="column" height="100%" width="100%">
    <Box flex="1 1 auto">
      <Plot
        data={finalData}
        layout={finalLayout}
        floatTitle="Trajectory"
        fixedMargins={true}
        className={styles.trajectory}
        data-testid='trajectory'
        metaInfoLink={archiveURL}
        {...other}
      />
    </Box>
    {methodology && <Box flex="0 0 auto">
      <PropertyMethodologyList xs={12}>
        <PropertyMethodologyItem
          title="Molecular dynamics"
          data={methodology.molecular_dynamics}
          path={([...trajectoryPath, 'methodology', 'molecular_dynamics']).join('.')}
        />
      </PropertyMethodologyList>
    </Box>}
  </Box>
})

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
  className: PropTypes.string,
  classes: PropTypes.object,
  archiveURL: PropTypes.string // Path for the data in the archive browser
}

export default withErrorHandler(trajectoryError)(Trajectory)

/**
 * Fetches all trajectories from the archive and displays them.
 */
const useTrajectoriesStyles = makeStyles((theme) => ({
  trajectory: {height: '450px'}
}))
const TrajectoriesRaw = React.memo(({index, archive}) => {
  const styles = useTrajectoriesStyles()
  const urlPrefix = `${getLocation()}/data`

  // Resolve and return component for showing the list of trajectories
  const trajsIndex = get(index, trajectoryPath) || []
  const trajsArchive = get(archive, trajectoryPath) || []

  return trajsIndex.map((trajIndex, i) => {
    const trajProperties = new Set(trajIndex?.available_properties || [])
    const trajArchive = trajsArchive[i]
    const hasPressure = trajProperties.has('pressure')
    const pressure = hasPressure ? trajArchive?.pressure : false
    const hasTemperature = trajProperties.has('temperature')
    const temperature = hasTemperature ? trajArchive?.temperature : false
    const hasEnergyPotential = trajProperties.has('energy_potential')
    const energyPotential = hasEnergyPotential ? trajArchive?.energy_potential : false
    const methodology = trajIndex?.methodology || false

    return <Trajectory
      key={i}
      pressure={pressure}
      temperature={temperature}
      energyPotential={energyPotential}
      methodology={methodology}
      classes={{trajectory: styles.trajectory}}
      archiveURL={`${urlPrefix}/${trajectoryPath.join('/')}:${i}`}
    />
  })
})

TrajectoriesRaw.propTypes = {
  index: PropTypes.object,
  archive: PropTypes.object
}

export const Trajectories = withErrorHandler(trajectoryError)(TrajectoriesRaw)
