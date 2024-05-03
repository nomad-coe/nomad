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
import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import { capitalize, isEmpty } from 'lodash'
import { Select, FormControl, MenuItem } from '@material-ui/core'
import { PropertyCard, PropertyGrid, PropertyItem, PropertyCardActions } from './PropertyCard'
import Quantity, { QuantityTable, QuantityRow, QuantityCell } from '../../Quantity'
import { MaterialButton } from '../../nav/Routes'
import Structure, { toMateriaStructure } from '../../visualization/Structure'
import NoData from '../../visualization/NoData'
import { encyclopediaBase, guiBase } from '../../../config'

/**
 * For displaying the most descriptive chemical formula that is present in an
 * entry.
*/
export function Formula({data}) {
  const formula = (data) => {
    const material = data?.results?.material
    if (!material) {
      return null
    }

    return material.chemical_formula_hill || material.chemical_formula_reduced || material.chemical_formula_descriptive
  }

  return <Quantity
    quantity={formula} label='Formula' noWrap data={data}
    description="The chemical formula that describes the simulated system or experimental sample."
  />
}

Formula.propTypes = {
  data: PropTypes.object
}

/**
 * Displays a summary of material related properties for an entry.
 */
const MaterialCard = React.memo(({index, properties, archive}) => {
  // Find out which properties are present
  const structures = index?.results?.properties?.structures
  const hasStructures = structures?.structure_original ||
    structures?.structure_conventional ||
    structures?.structure_primitive
  const hasLatticeParameters = structures?.structure_original?.lattice_parameters ||
    structures?.structure_conventional?.lattice_parameters ||
    structures?.structure_primitive?.lattice_parameters
  const hasSymmetry = index?.results?.material?.symmetry

  // Get all structure types that are present
  const structureTypes = {}
  if (structures) {
    Object.keys(structures).forEach((key) => {
      const name = capitalize(key.split('_').pop())
      structureTypes[key] = name
    })
  }

  // Determine the structure to show by default
  let defaultStructure
  if (structures && !isEmpty(structures)) {
    if ('structure_original' in structures) {
      defaultStructure = 'structure_original'
    } else if ('structure_conventional' in structures) {
      defaultStructure = 'structure_conventional'
    } else if (!isEmpty(structures)) {
      defaultStructure = Object.keys(structures)[0]
    }
  }
  const [structureType, setStructureType] = useState(defaultStructure)

  const handleStructureChange = useCallback((event) => {
    setStructureType(event.target.value)
  }, [])

  if (!index?.results?.material) {
    return null
  }

  // Prepare the data for the visualizer
  const urlPrefix = `${window.location.pathname.slice(guiBase.length)}/data/results/properties/structures`
  const materialId = index.results?.material?.material_id
  const structurePath = `results.properties.structures.${structureType}`
  const structureSection = archive?.results?.properties?.structures?.[structureType]
  const m_path = `${urlPrefix}/${structureType}`
  const structure = structureSection && toMateriaStructure(structureSection)

  // Dropdown for selecting a specific structure
  const select = hasStructures && <FormControl>
    <Select
      value={structureType}
      onChange={handleStructureChange}
      label="Structure"
    >
      {Object.entries(structureTypes).map(([type, name]) =>
        <MenuItem key={type} value={type}>{name}</MenuItem>
      )}
    </Select>
  </FormControl>
  return <PropertyCard title="Material" action={select}>
    <PropertyGrid>
      <PropertyItem title="Composition" xs={6} height="auto">
        <QuantityTable data={index}>
          <QuantityRow>
            <QuantityCell colSpan={2}>
              <Formula data={index}/>
            </QuantityCell>
          </QuantityRow>
          <QuantityRow>
            <QuantityCell quantity="results.material.structural_type"/>
          </QuantityRow>
          <QuantityRow>
            <QuantityCell quantity="results.material.elements" colSpan={2}/>
          </QuantityRow>
          <QuantityRow>
            <QuantityCell quantity="results.material.n_elements"/>
          </QuantityRow>
        </QuantityTable>
      </PropertyItem>
      <PropertyItem title="Structure" xs={6} height="280px">
        {hasStructures
          ? <Structure
            data={structure}
            structuralType={index.results?.material?.structural_type}
            cellType={structureType.split('_').pop()}
            m_path={m_path}
            data-testid="viewer-material"
          />
          : <NoData/>}
      </PropertyItem>
      <PropertyItem title="Symmetry" xs={6} height="auto" minHeight="200px">
        {hasSymmetry
          ? <QuantityTable data={index}>
            <QuantityRow>
              <QuantityCell quantity="results.material.symmetry.crystal_system"/>
              <QuantityCell quantity="results.material.symmetry.bravais_lattice"/>
            </QuantityRow>
            <QuantityRow>
              <QuantityCell quantity="results.material.symmetry.space_group_number"/>
              <QuantityCell quantity="results.material.symmetry.space_group_symbol"/>
            </QuantityRow>
            <QuantityRow>
              <QuantityCell quantity="results.material.symmetry.point_group"/>
              <QuantityCell quantity="results.material.symmetry.structure_name"/>
            </QuantityRow>
          </QuantityTable>
          : <NoData/>
        }
      </PropertyItem>
      <PropertyItem title="Lattice parameters" xs={6} height="auto" minHeight="200px">
        {hasLatticeParameters
          ? <QuantityTable data={index}>
            <QuantityRow>
              <QuantityCell quantity={`${structurePath}.lattice_parameters.a`}/>
              <QuantityCell quantity={`${structurePath}.lattice_parameters.b`}/>
              <QuantityCell quantity={`${structurePath}.lattice_parameters.c`}/>
            </QuantityRow>
            <QuantityRow>
              <QuantityCell quantity={`${structurePath}.lattice_parameters.alpha`} label="α"/>
              <QuantityCell quantity={`${structurePath}.lattice_parameters.beta`} label="β"/>
              <QuantityCell quantity={`${structurePath}.lattice_parameters.gamma`} label="γ"/>
            </QuantityRow>
            <QuantityRow>
              <QuantityCell colSpan={2} quantity={`${structurePath}.cell_volume`}/>
            </QuantityRow>
          </QuantityTable>
          : <NoData/>}
      </PropertyItem>
    </PropertyGrid>
    {encyclopediaBase && materialId &&
      <PropertyCardActions>
        <MaterialButton materialId={materialId}>
          View in Encyclopedia
        </MaterialButton>
      </PropertyCardActions>
    }
  </PropertyCard>
})

MaterialCard.propTypes = {
  index: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default MaterialCard
