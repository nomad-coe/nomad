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
import { useUnits } from '../../../units'
import { resolveRef } from '../../archive/metainfo'
import GeometryOptimization from '../../visualization/GeometryOptimization'

export default function GeometryOptimizationCard({index, archive, properties}) {
  const units = useUnits()

  // Find out which properties are present
  const hasGeometryOptimization = properties.has('geometry_optimization')

  // Do not show the card if none of the properties are available, or if only
  // one step is calculated.
  if (!hasGeometryOptimization) {
    return null
  }

  // Resolve energies
  let energies = hasGeometryOptimization ? null : false
  const geoOptPropsArchive = archive?.results?.properties?.geometry_optimization
  if (hasGeometryOptimization && archive) {
    energies = resolveRef(geoOptPropsArchive.energies, archive)
  }

  // Resolve convergence properties
  let convergence = false
  const geoOptProps = index?.results?.properties?.geometry_optimization
  const geoOptMethod = index.results?.method?.simulation?.geometry_optimization
  if (hasGeometryOptimization) {
    convergence = {...geoOptMethod, ...geoOptProps}
  }

  return <PropertyCard title="Geometry optimization">
    <GeometryOptimization
      energies={energies}
      convergence={convergence}
      units={units}
    />
  </PropertyCard>
}

GeometryOptimizationCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}
