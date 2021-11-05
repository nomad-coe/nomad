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
import EELS from '../../visualization/EELS'
import { resolveRef } from '../../archive/metainfo'

/**
 * Card displaying spectroscopic properties.
*/
const SpectroscopyCard = React.memo(({index, properties, archive}) => {
  const units = useUnits()

  // Find out which properties are present
  const hasEELS = properties.has('eels')

  // Do not show the card if none of the properties are available
  if (!hasEELS) {
    return null
  }

  // Resolve EELS data
  const dataRef = archive?.results?.properties?.spectroscopy?.eels
  const data = dataRef && resolveRef(dataRef, archive)

  return <PropertyCard title="Spectroscopy">
    <PropertyGrid>
      <PropertyItem title="Electron energy loss spectroscopy" xs={12}>
        <EELS
          data={data}
          layout={{yaxis: {autorange: true}}}
          aspectRatio={2}
          units={units}
        />
      </PropertyItem>
    </PropertyGrid>
  </PropertyCard>
})

SpectroscopyCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default SpectroscopyCard
