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
import React from 'react'
import PropTypes from 'prop-types'
import { PropertyCard } from './PropertyCard'
import { resolveInternalRef } from '../../../utils'
import GeometryOptimization from '../../visualization/GeometryOptimization'

const GeometryOptimizationCard = React.memo(({index, archive, properties}) => {
  const geoOptProps = index?.results?.properties?.geometry_optimization

  // Find out which properties are present. If only one step is calculated
  // (n_calculations == 1 and none of the convergence criteria are available),
  // the card will not be displayed.
  const hasEnergies = properties.has('geometry_optimization') &&
    index?.results?.properties?.n_calculations > 1
  const hasConvergence = properties.has('geometry_optimization') &&
    (
      index?.results?.properties?.geometry_optimization?.final_energy_difference ||
      index?.results?.properties?.geometry_optimization?.final_displacement_maximum ||
      index?.results?.properties?.geometry_optimization?.final_force_maximum
    )

  // Do not show the card if none of the properties are available
  if (!hasEnergies && !hasConvergence) {
    return null
  }

  // Resolve energies
  let energies = hasEnergies ? null : false
  const energiesArchive = archive?.results?.properties?.geometry_optimization?.energies
  if (hasEnergies && energiesArchive) {
    energies = resolveInternalRef(energiesArchive, archive)
  }

  // Resolve convergence properties
  const convergence = hasConvergence ? geoOptProps : false

  return <PropertyCard title="Geometry optimization">
    <GeometryOptimization
      energies={energies}
      convergence={convergence}
    />
  </PropertyCard>
})

GeometryOptimizationCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default GeometryOptimizationCard
