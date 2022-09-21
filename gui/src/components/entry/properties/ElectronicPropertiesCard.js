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
import { getLocation, resolveInternalRef } from '../../../utils'
import ElectronicProperties from '../../visualization/ElectronicProperties'
import { refPath } from '../../archive/metainfo'

const ElectronicPropertiesCard = React.memo(({index, properties, archive}) => {
  const units = useUnits()
  const urlPrefix = `${getLocation()}/data`

  // Find out which properties are present
  const hasDos = properties.has('dos_electronic')
  const hasBs = properties.has('band_structure_electronic')

  // Do not show the card if none of the properties are available
  if (!hasDos && !hasBs) {
    return null
  }

  // Resolve DOS data
  let dos = hasDos ? null : false
  const dosData = archive?.results?.properties?.electronic?.dos_electronic
  if (dosData) {
    dos = {}
    dos.energies = resolveInternalRef(dosData.energies, archive)
    dos.densities = resolveInternalRef(dosData.total, archive).map(dos => dos.value)
    if (dosData.band_gap) {
      dos.energy_highest_occupied = Math.max(...dosData.band_gap.map(x => x.energy_highest_occupied))
    }
    dos.m_path = `${urlPrefix}/${refPath(dosData.energies.split('/').slice(0, -1).join('/'))}`
  }

  // Resolve band structure data
  let bs = hasBs ? null : false
  const bsData = archive?.results?.properties?.electronic?.band_structure_electronic
  if (bsData) {
    bs = {}
    bs.reciprocal_cell = resolveInternalRef(bsData.reciprocal_cell, archive)
    bs.segment = resolveInternalRef(bsData.segment, archive)
    if (bsData.band_gap) {
      bs.energy_highest_occupied = Math.max(...bsData.band_gap.map(x => x.energy_highest_occupied))
      bs.band_gap = bsData.band_gap
    }
    bs.m_path = `${urlPrefix}/${refPath(bsData.reciprocal_cell.split('/').slice(0, -1).join('/'))}`
  }

  return <PropertyCard title="Electronic properties">
    <ElectronicProperties bs={bs} dos={dos} units={units} />
  </PropertyCard>
})

ElectronicPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default ElectronicPropertiesCard
