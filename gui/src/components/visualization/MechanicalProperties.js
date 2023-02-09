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
import React, { } from 'react'
import PropTypes from 'prop-types'
import { SectionTable } from '../Quantity'
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'
import EnergyVolumeCurve from '../visualization/EnergyVolumeCurve'

/**
 * Shows a summary of Mechanical properties.
 */
const modulusQuantities = {
  type: {label: 'Type', align: 'left'},
  value: {label: 'Value', align: 'right'}
}
const MechanicalProperties = React.memo(({
  evCurves,
  bulkModulus,
  shearModulus,
  units
}) => {
  return <PropertyGrid>
    {evCurves !== false && <PropertyItem title="Energy-volume curve" xs={12} height="500px">
      <EnergyVolumeCurve
        data={evCurves}
        units={units}
        data-testid="energy-volume-curve"
      />
    </PropertyItem>}
    <PropertyItem title="Bulk modulus" xs={6} height="auto">
      <SectionTable
        horizontal
        section="results.properties.mechanical.bulk_modulus"
        quantities={modulusQuantities}
        data={bulkModulus}
        units={units}
        data-testid="bulk-modulus"
      />
    </PropertyItem>
    <PropertyItem title="Shear modulus" xs={6} height="auto">
      <SectionTable
        horizontal
        section="results.properties.mechanical.shear_modulus"
        quantities={modulusQuantities}
        data={shearModulus}
        units={units}
        data-testid="shear-modulus"
      />
    </PropertyItem>
  </PropertyGrid>
})

MechanicalProperties.propTypes = {
  evCurves: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  bulkModulus: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  shearModulus: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  units: PropTypes.object // Contains the unit configuration
}

export default MechanicalProperties
