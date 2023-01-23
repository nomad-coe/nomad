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
import StructuralProperties from '../../visualization/StructuralProperties'

/**
 * Card for displaying structural properties.
 */
const StructuralPropertiesCard = React.memo(({index, properties, archive}) => {
  // Check what data is available and do not show the card if none of the properties are
  // available.
   // Find out which properties are present
  const hasRdf = properties.has('radial_distribution_function')
  const hasRg = properties.has('radius_of_gyration')

  // Do not show the card if none of the properties are available
  if (!hasRdf && !hasRg) {
    return null
  }

  return <PropertyCard title="Structural properties">
    <StructuralProperties
      index={index}
      properties={properties}
      archive={archive}
    />
  </PropertyCard>
})

StructuralPropertiesCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default StructuralPropertiesCard
