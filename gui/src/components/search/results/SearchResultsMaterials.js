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
import React, { useCallback } from 'react'
import PropTypes from 'prop-types'
import { IconButton } from '@material-ui/core'
import DetailsIcon from '@material-ui/icons/MoreHoriz'
import searchQuantities from '../../../searchQuantities'
import { encyclopediaEnabled } from '../../../config'
import { MaterialButton } from '../../nav/Routes'
import NewDataTable from '../../NewDataTable'

/**
 * Displays the list of search results for materials.
 */
const columns = {
  formula: {
    label: 'Formula',
    render: row => row.chemical_formula_hill,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.chemical_formula_hill'].description
  },
  structural_type: {
    label: 'Structural type',
    render: row => row.structural_type,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.structural_type'].description
  },
  structure_name: {
    label: 'Structure name',
    render: row => row.symmetry?.structure_name,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.structure_name'].description
  },
  crystal_system: {
    label: 'Crystal system',
    render: row => row.symmetry?.crystal_system,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.crystal_system'].description
  },
  space_group_number: {
    label: 'Space group number',
    render: row => row.symmetry?.space_group_number,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.space_group_number'].description
  },
  space_group_symbol: {
    label: 'Space group symbol',
    render: row => row.symmetry?.space_group_symbol,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.symmetry.space_group_symbol'].description
  },
  material_id: {
    label: 'Material ID',
    render: row => row.material_id,
    supportsSort: true,
    ellipsisFront: true,
    description: searchQuantities['results.material.material_id'].description
  }
}
const selectedColumns = ['formula', 'structural_type', 'structure_name', 'space_group_number', 'crystal_system', 'material_id']

const SearchResultsMaterials = React.memo(({
  data,
  query,
  className,
  ...rest
}) => {
  const total = data?.pagination && data.pagination.total

  const renderActions = useCallback((row, selected) => (
    <MaterialButton
      materialId={row.material_id}
      component={IconButton}
    >
      <DetailsIcon/>
    </MaterialButton>
  ), [])

  return <NewDataTable
    entityLabels={['material', 'materials']}
    id={row => row.material_id}
    total={total}
    columns={columns}
    selectedColumns={selectedColumns}
    selectedColumnsKey="materials"
    entryActions={encyclopediaEnabled ? renderActions : undefined}
    data={data?.data || []}
    rows={data?.data.length || 0}
    {...rest}
  />
})
SearchResultsMaterials.propTypes = {
  data: PropTypes.object,
  query: PropTypes.object,
  className: PropTypes.string
}

export default SearchResultsMaterials
