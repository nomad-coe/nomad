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
import { QuantityTable, QuantityRow, QuantityCell } from '../../Quantity'

/**
 * Displays a summary of material related properties for an entry.
 */
const MaterialCardFormula = React.memo(({index}) => {
  if (!index?.results?.material) {
    return null
  }

  return <PropertyCard title="Material">
    <PropertyGrid>
      <PropertyItem title="Composition" xs={12} height="auto">
        <QuantityTable data={index}>
          <QuantityRow>
            <QuantityCell quantity="results.material.chemical_formula_iupac"/>
            <QuantityCell quantity="results.material.chemical_formula_reduced"/>
          </QuantityRow>
          <QuantityRow>
            <QuantityCell quantity="results.material.chemical_formula_hill"/>
            <QuantityCell quantity="results.material.chemical_formula_descriptive"/>
          </QuantityRow>
          <QuantityRow>
            <QuantityCell quantity="results.material.elements"/>
            <QuantityCell quantity="results.material.n_elements"/>
          </QuantityRow>
        </QuantityTable>
      </PropertyItem>
    </PropertyGrid>
  </PropertyCard>
})

MaterialCardFormula.propTypes = {
  index: PropTypes.object
}

export default MaterialCardFormula
