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
import { resolveInternalRef } from '../../../utils'
import { PropertyCard } from './PropertyCard'
import ElectronicProperties from '../../visualization/ElectronicProperties'
import { resolveGreensFunctions } from '../../visualization/GreensFunctions'
import { resolveDosNew, resolveDosOld } from '../../visualization/DOS'

const ElectronicPropertiesCard = React.memo(({index, properties, archive}) => {
  // Find out which properties are present
  const hasDosNew = properties.has('dos_electronic_new')
  const hasDosOld = properties.has('dos_electronic')
  const hasBs = properties.has('band_structure_electronic')
  const hasBandGap = properties.has('electronic.band_structure_electronic.band_gap')
  const hasGf = properties.has('greens_functions_electronic')

  // Do not show the card if none of the properties are available
  if (!hasDosNew && !hasDosOld && !hasBs && !hasBandGap && !hasGf) return null

  let bsReferences = archive?.results?.properties?.electronic?.band_structure_electronic || []
  const pattern = '\\.\\./(?:entries|upload/archive|uploads.+?archive)/(.+?)(?:/archive#|#)(.+)'
  if (!Array.isArray(bsReferences)) bsReferences = [bsReferences]

  // Resolve the new DOS schema data. For older entries, resolve the old DOS schema data
  let dos = (hasDosOld || hasDosNew) ? undefined : false
  if (archive) {
    if (hasDosNew) dos = resolveDosNew(properties, archive, pattern)
    if (!hasDosNew && hasDosOld) dos = resolveDosOld(properties, archive, pattern)
  }

  // Resolve band structure data
  let bs = hasBs ? undefined : false
  let bz = hasBs ? undefined : false
  if (archive) {
    bs = []
    try {
      for (const reference of bsReferences) {
        const d = {}
        const match = reference.reciprocal_cell?.match(pattern)
        const path = match ? match[2] : reference.reciprocal_cell
        const segmentPath = match ? reference.segment.map(ref => ref.match(pattern)[2]) : reference.segment
        const sourceArchive = match ? (archive.m_ref_archives[match[1]] || archive.m_ref_archives[reference.reciprocal_cell.split('#')[0]]) : archive
        d.segment = sourceArchive ? resolveInternalRef(segmentPath, sourceArchive) : null
        d.name = reference.label
        d.m_path = `${archive?.metadata?.entry_id}/data/results/properties/electronic/band_structure_electronic`
        if (reference.band_gap) {
          d.energy_highest_occupied = Math.max(...reference.band_gap.map(x => x.energy_highest_occupied))
          d.band_gap = reference.band_gap
        }
        if (d.segment) bs.push(d)
        const reciprocal_cell = sourceArchive ? resolveInternalRef(path, sourceArchive) : null
        if (reciprocal_cell) {
          bz = {
            reciprocal_cell: reciprocal_cell,
            segment: d.segment
          }
        }
      }
    } catch (e) {
    }
    bs = bs.length === 0 ? false : bs
    bz = bz || false
  }

  // Resolve band gap data. TODO: The API is not returning band gap
  // information when using electronic: 'include-resolved' and there is
  // nothing to resolve. This is why we alternatively also look for the band
  // gap data in the index data.
  let bg
  if (hasBandGap) {
    bg = []
    function addBandGaps(sections) {
      for (const section of sections || []) {
        if (!section.band_gap) continue
        bg.push(...section.band_gap)
      }
    }
    addBandGaps(bsReferences)
    if (bg.length === 0) {
      let bsReferencesIndex = index?.results?.properties?.electronic?.band_structure_electronic || []
      if (!Array.isArray(bsReferencesIndex)) bsReferencesIndex = [bsReferencesIndex]
      addBandGaps(bsReferencesIndex)
    }
    bg = bg.length === 0 ? false : bg
  } else {
    bg = false
  }

  // Resolve Greens functions data
  const gf = resolveGreensFunctions(properties, archive, pattern)

  return <PropertyCard title="Electronic properties">
    <ElectronicProperties
      bs={bs}
      dos={dos}
      brillouin_zone={bz}
      band_gap={bg}
      gf={gf}
      index={index}
    />
  </PropertyCard>
})

ElectronicPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default ElectronicPropertiesCard
