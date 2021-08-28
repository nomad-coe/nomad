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
import PropertyCard from './PropertyCard'
import { useUnits } from '../../../units'
import ElectronicProperties from '../../visualization/ElectronicProperties'
import { refPath, resolveRef } from '../../archive/metainfo'

export default function ElectronicPropertiesCard({entryMetadata, archive}) {
  const units = useUnits()

  const properties = new Set(
    entryMetadata.results ? entryMetadata.results.properties.available_properties : [])
  const hasDos = properties.has('dos_electronic')
  const hasBs = properties.has('band_structure_electronic')

  if (!hasDos && !hasBs) {
    return null
  }

  const archiveUrl = `/entry/id/${entryMetadata.upload_id}/${entryMetadata.entry_id}/archive`

  let dos = hasDos ? null : false
  const dosData = archive?.results?.properties?.electronic?.dos_electronic
  if (dosData) {
    dos = {}
    dos.energies = resolveRef(dosData.energies, archive)
    dos.densities = resolveRef(dosData.densities, archive)
    if (dosData.channel_info) {
      dos.energy_highest_occupied = Math.max(...dosData.channel_info.map(x => x.energy_highest_occupied))
    }
    dos.m_path = `${archiveUrl}/${refPath(dosData.energies.split('/').slice(0, -1).join('/'))}`
  }

  let bs = hasBs ? null : false
  const bsData = archive?.results?.properties?.electronic?.band_structure_electronic
  if (bsData) {
    bs = {}
    bs.reciprocal_cell = resolveRef(bsData.reciprocal_cell, archive)
    bs.segments = resolveRef(bsData.segments, archive)
    if (bsData.channel_info) {
      bs.energy_highest_occupied = Math.max(...bsData.channel_info.map(x => x.energy_highest_occupied))
      bs.channel_info = bsData.channel_info
    }
    bs.m_path = `${archiveUrl}/${refPath(bsData.reciprocal_cell.split('/').slice(0, -1).join('/'))}`
  }

  return <PropertyCard title="Electronic properties">
    <ElectronicProperties bs={bs} dos={dos} units={units} />
  </PropertyCard>
}

ElectronicPropertiesCard.propTypes = {
  entryMetadata: PropTypes.object.isRequired,
  archive: PropTypes.object
}
