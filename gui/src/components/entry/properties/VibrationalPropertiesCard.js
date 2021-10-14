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
import { refPath, resolveRef } from '../../archive/metainfo'
import VibrationalProperties from '../../visualization/VibrationalProperties'

/**
 * Card displaying vibrational properties.
 */
const VibrationalPropertiesCard = React.memo(({entryMetadata, properties, archive}) => {
  const units = useUnits()

  // Find out which properties are present
  const hasDos = properties.has('dos_phonon')
  const hasBs = properties.has('band_structure_phonon')
  const hasEnergyFree = properties.has('energy_free_helmholtz')
  const hasHeatCapacity = properties.has('heat_capacity_constant_volume')

  // Do not show the card if none of the properties are available
  if (!hasDos && !hasBs && !hasEnergyFree && !hasHeatCapacity) {
    return null
  }

  const archiveUrl = `/entry/id/${entryMetadata.upload_id}/${entryMetadata.entry_id}/archive`

  let dos = hasDos ? null : false
  const dosData = archive?.results?.properties?.vibrational?.dos_phonon
  if (dosData) {
    dos = {}
    dos.energies = resolveRef(dosData.energies, archive)
    dos.densities = resolveRef(dosData.total, archive).map(dos => dos.value)
    dos.m_path = `${archiveUrl}/${refPath(dosData.energies.split('/').slice(0, -1).join('/'))}`
  }

  let bs = hasBs ? null : false
  const bsData = archive?.results?.properties?.vibrational?.band_structure_phonon
  if (bsData) {
    bs = {}
    bs.segment = resolveRef(bsData.segment, archive)
    bs.m_path = `${archiveUrl}/${refPath(bsData.segment[0].split('/').slice(0, -2).join('/'))}`
  }

  let energyFree = hasEnergyFree ? null : false
  const energyFreeData = archive?.results?.properties?.vibrational?.energy_free_helmholtz
  if (energyFreeData) {
    energyFree = {}
    energyFree.energies = resolveRef(energyFreeData.energies, archive)
    energyFree.temperatures = resolveRef(energyFreeData.temperatures, archive)
    energyFree.m_path = `${archiveUrl}/${refPath(energyFreeData.energies)}`
  }

  let heatCapacity = hasHeatCapacity ? null : false
  const heatCapacityData = archive?.results?.properties?.vibrational?.heat_capacity_constant_volume
  if (heatCapacityData) {
    heatCapacity = {}
    heatCapacity.heat_capacities = resolveRef(heatCapacityData.heat_capacities, archive)
    heatCapacity.temperatures = resolveRef(heatCapacityData.temperatures, archive)
    heatCapacity.m_path = `${archiveUrl}/${refPath(heatCapacityData.heat_capacities)}`
  }

  return <PropertyCard title="Vibrational properties">
    <VibrationalProperties
      bs={bs}
      dos={dos}
      heatCapacity={heatCapacity}
      freeEnergy={energyFree}
      units={units}
    />
  </PropertyCard>
})

VibrationalPropertiesCard.propTypes = {
  entryMetadata: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default VibrationalPropertiesCard
