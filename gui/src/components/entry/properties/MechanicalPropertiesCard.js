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
import { useUnitContext } from '../../units/UnitContext'
import { getLocation, resolveInternalRef } from '../../../utils'
import { refPath } from '../../archive/metainfo'
import MechanicalProperties from '../../visualization/MechanicalProperties'

/**
 * Card displaying mechanical properties.
 */
const MechanicalPropertiesCard = React.memo(({index, properties, archive}) => {
  const {units} = useUnitContext()
  const urlPrefix = `${getLocation()}/data`

  // Find out which properties are present
  const hasEVCurves = properties.has('energy_volume_curve')
  const hasBulkModulus = properties.has('bulk_modulus')
  const hasShearModulus = properties.has('shear_modulus')

  // Do not show the card if none of the properties are available
  if (!hasBulkModulus && !hasEVCurves && !hasShearModulus) {
    return null
  }

  // Resolve EV curves
  let evCurves = hasEVCurves ? null : false
  if (hasEVCurves && archive) {
    const evCurveData = archive?.results?.properties?.mechanical?.energy_volume_curve
    evCurves = {
      m_path: `${urlPrefix}/${refPath(evCurveData[0].volumes.split('/').slice(0, -1).join('/'))}`,
      data: evCurveData.map(ev => ({
        volumes: resolveInternalRef(ev.volumes, archive),
        energies: resolveInternalRef(ev.energies_raw || ev.energies_fit, archive),
        name: ev.type
      }))
    }
  }

  // Resolve bulk modulus
  let bulkModulus = false
  if (hasBulkModulus) {
    bulkModulus = {
      data: index?.results?.properties?.mechanical?.bulk_modulus
    }
  }

  // Resolve shear modulus
  let shearModulus = false
  if (hasShearModulus) {
    shearModulus = {
      data: index?.results?.properties?.mechanical?.shear_modulus
    }
  }

  return <PropertyCard title="Mechanical properties">
    <MechanicalProperties
      evCurves={evCurves}
      bulkModulus={bulkModulus}
      shearModulus={shearModulus}
      units={units}
    />
  </PropertyCard>
})

MechanicalPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default MechanicalPropertiesCard
