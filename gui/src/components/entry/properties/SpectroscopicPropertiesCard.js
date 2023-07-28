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
import { PropertyCard, PropertyGrid, PropertyItem } from './PropertyCard'
import Spectra from '../../visualization/Spectra'

/**
 * Card displaying spectroscopic properties.
*/
const SpectroscopicPropertiesCard = React.memo(({index, properties, archive}) => {
  // Find out which properties are present
  const hasSpectra = properties.has('spectra')

  // Do not show the card if none of the properties are available
  if (!hasSpectra) return null

  // Check index for which types of spectras exist. Create a mapping based on
  // this: this will allow the GUI to show a placeholder while the data is
  // loaded from archive.
  const spectras = hasSpectra ? {} : false
  for (const spectra of index?.results?.properties?.spectroscopic?.spectra) {
    spectras[spectra.type] = undefined
  }

  // When archive is loaded, fill out the information
  if (archive) {
    let spectraReferences = archive?.results?.properties?.spectroscopic?.spectra
    if (!Array.isArray(spectraReferences)) spectraReferences = [spectraReferences]
    if (spectraReferences) {
      for (const spectra of spectraReferences) {
        const data = {}
        data.type = spectra.type
        data.label = spectra.label === 'computation' ? 'comp.' : 'exp.'
        data.intensities = spectra.intensities
        data.energies = spectra.energies
        if (data.energies) {
          if (spectras[spectra.type]) {
            spectras[spectra.type].push(data)
          } else {
            spectras[spectra.type] = [data]
          }
        }
      }
    }
  }

  return <PropertyCard title="Spectroscopic properties">
    <PropertyGrid>
    {spectras && Object.entries(spectras).map(([key, value]) => {
      return <PropertyItem key={key} title={key} xs={12}>
        <Spectra
          data={value}
          layout={{yaxis: {autorange: true}}}
        />
      </PropertyItem>
    })}
    </PropertyGrid>
  </PropertyCard>
})

SpectroscopicPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default SpectroscopicPropertiesCard
