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
import React, { useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { resolveInternalRef } from '../../../utils'
import { PropertyCard } from './PropertyCard'
import ElectronicProperties from '../../visualization/ElectronicProperties'

const ElectronicPropertiesCard = React.memo(({index, properties, archive}) => {
  const [data, setData] = useState([false, false, false, false])
  const [refArchives, setRefArchives] = useState()

  // Find out which properties are present
  const hasDos = properties.has('dos_electronic')
  const hasBs = properties.has('band_structure_electronic')

  // TODO remove refArchives after fixing archive.m_ref_archives not vanishing when switching to data tab
  useEffect(() => {
    if (Object.keys(archive.m_ref_archives).length > 0) {
      const archives = {}
      // TODO remove this when the path format for the entry is fixed
      // it changes from ../uploads/upload_id/archive/entry_id to ../uploa/archive/entry_id
      const pattern = '\\.\\./(?:entries|upload/archive|uploads.+?archive)/(.+)'
      Object.keys(archive.m_ref_archives).forEach(path => {
        const match = path.match(pattern)
        if (match) archives[match[1]] = archive.m_ref_archives[path]
      })
      setRefArchives(archives)
    }
  }, [archive])

  archive.m_ref_archives = refArchives || archive.m_ref_archives

  // TODO remove useEffect once archive.m_ref_archives is fixed
  useEffect(() => {
    const dosReferences = archive?.results?.properties?.electronic?.dos_electronic || []
    const bsReferences = archive?.results?.properties?.electronic?.band_structure_electronic || []
    const pattern = '\\.\\./(?:entries|upload/archive|uploads.+?archive)/(.+?)(?:/archive#|#)(.+)'

    // Resolve DOS data
    let references = dosReferences
    if (!Array.isArray(references)) references = [references]
    let dos = []
    for (const reference of references) {
      const d = {}
      const match = reference.energies.match(pattern)
      const path = match ? match[2] : reference.energies
      const totalPath = match ? reference.total.map(ref => ref.match(pattern)[2]) : reference.total
      const sourceArchive = match ? (archive.m_ref_archives[match[1]] || archive.m_ref_archives[reference.energies.split('#')[0]]) : archive
      if (sourceArchive) {
        d.energies = resolveInternalRef(path, sourceArchive)
        const internalRef = resolveInternalRef(totalPath, sourceArchive)
        d.densities = internalRef.map(dos => dos.value)
        d.normalization_factors = internalRef.map(dos => dos.normalization_factor)
      }
      d.name = reference.label
      if (reference.band_gap) {
        d.energy_highest_occupied = Math.max(...reference.band_gap.map(x => x.energy_highest_occupied))
      }
      d.m_path = `${archive?.metadata?.entry_id}/data/results/properties/electronic/dos_electronic`
      if (d.energies && d.densities) dos.push(d)
    }
    dos = dos.length === 0 ? false : dos

    // Resolve band structure data
    references = bsReferences
    if (!Array.isArray(references)) references = [references]
    let bs = []
    let band_gap = []
    let brillouin_zone = false
    brillouin_zone = false
    for (const reference of references) {
      const d = {}
      const match = reference.reciprocal_cell.match(pattern)
      const path = match ? match[2] : reference.reciprocal_cell
      const segmentPath = match ? reference.segment.map(ref => ref.match(pattern)[2]) : reference.segment
      const sourceArchive = match ? (archive.m_ref_archives[match[1]] || archive.m_ref_archives[reference.reciprocal_cell.split('#')[0]]) : archive
      d.segment = sourceArchive ? resolveInternalRef(segmentPath, sourceArchive) : null
      d.name = reference.label
      d.m_path = `${archive?.metadata?.entry_id}/data/results/properties/electronic/band_structure_electronic`
      if (reference.band_gap) {
        d.energy_highest_occupied = Math.max(...reference.band_gap.map(x => x.energy_highest_occupied))
        d.band_gap = reference.band_gap
        band_gap.push(...reference.band_gap)
      }
      if (d.segment) bs.push(d)
      const reciprocal_cell = sourceArchive ? resolveInternalRef(path, sourceArchive) : null
      if (reciprocal_cell) {
        brillouin_zone = {
          reciprocal_cell: reciprocal_cell,
          segment: d.segment
        }
      }
    }

    bs = bs.length === 0 ? false : bs
    band_gap = band_gap.length === 0 ? false : band_gap
    setData([dos, bs, brillouin_zone, band_gap])
  }, [archive])

  // Do not show the card if none of the properties are available
  if (!hasDos && !hasBs) {
    return null
  }

  return <PropertyCard title="Electronic properties">
    <ElectronicProperties
      bs={data[1]}
      dos={data[0]}
      brillouin_zone={data[2]}
      band_gap={data[3]}
    />
  </PropertyCard>
})

ElectronicPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default ElectronicPropertiesCard
