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
import React, { useMemo } from 'react'
import PropTypes from 'prop-types'
import { PropertyCard } from './PropertyCard'
import { useUnits } from '../../../units'
import { getLocation } from '../../../utils'
import Trajectory from '../../visualization/Trajectory'

/**
 * Card displaying molecular dynamics properties.
 */
const ThermodynamicPropertiesCard = React.memo(({index, properties, archive}) => {
  const units = useUnits()
  const hasTrajectory = properties.has('trajectory')
  const urlPrefix = `${getLocation()}/data`

  // Do not show the card if none of the properties are available
  if (!hasTrajectory) {
    return null
  }

  // Resolve and return component for showing the list of trajectories
  const trajectory = useMemo(() => {
    const trajsIndex = index?.results?.properties?.thermodynamic?.trajectory
    const trajsArchive = archive?.results?.properties?.thermodynamic?.trajectory || []

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
        units={units}
        archiveURL={`${urlPrefix}/results/properties/thermodynamic/trajectory:${i}`}
      />
    })
  }, [archive, index, units, urlPrefix])

  return <PropertyCard title="Thermodynamic">
    {trajectory}
  </PropertyCard>
})

ThermodynamicPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default ThermodynamicPropertiesCard
