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
import { useUnits } from '../../../units'
import { SectionTable } from '../../Quantity'
import EELS from '../../visualization/EELS'
import { resolveInternalRef } from '../../../utils'

/**
 * Card displaying EELS properties.
*/
const eelsProperties = {
  detector_type: {label: 'Detector type', align: 'left'},
  resolution: {label: 'Resolution', align: 'right'},
  min_energy: {label: 'Min. energy', align: 'right'},
  max_energy: {label: 'Max. energy', align: 'right'}
}
const EELSPropertiesCard = React.memo(({index, properties, archive}) => {
  const units = useUnits()

  // Find out which properties are present
  const hasEELS = properties.has('eels')

  // Do not show the card if none of the properties are available
  if (!hasEELS) {
    return null
  }

  // Resolve EELS data
  let eelsTable
  let spectrumCurves
  if (archive) {
    const spectroscopy = archive?.results?.properties?.spectroscopy
    eelsTable = {data: [spectroscopy?.eels]}
    if (spectroscopy.spectrum) {
      spectrumCurves = [resolveInternalRef(spectroscopy.spectrum, archive)]
    }
  }

  return <PropertyCard title="EELS Properties">
    <PropertyGrid>
      <PropertyItem xs={12} height="25rem">
        <EELS
          data={spectrumCurves}
          layout={{yaxis: {autorange: true}}}
          units={units}
        />
      </PropertyItem>
      <PropertyItem xs={12} height="auto">
        <SectionTable
          horizontal
          section="results.properties.spectroscopy.eels"
          quantities={eelsProperties}
          data={eelsTable}
          units={units}
          data-testid="bulk-modulus"
          showIndex={spectrumCurves && spectrumCurves.length > 1}
        />
      </PropertyItem>
    </PropertyGrid>
  </PropertyCard>
})

EELSPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default EELSPropertiesCard
