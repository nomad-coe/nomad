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
import React, { useContext } from 'react'
import PropTypes from 'prop-types'
import { FilterSubMenu, filterMenuContext } from './FilterMenu'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import InputField from '../input/InputField'

const FilterSubMenuStructure = React.memo(({
  id,
  ...rest
}) => {
  const {selected, open} = useContext(filterMenuContext)
  const visible = open && id === selected

  return <FilterSubMenu id={id} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.structural_type"
          visible={visible}
          disableSearch
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.bravais_lattice"
          visible={visible}
          xs={6}
          disableSearch
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.crystal_system"
          visible={visible}
          xs={12}
          disableSearch
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.space_group_symbol"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.structure_name"
          visible={visible}
          initialSize={5}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.strukturbericht_designation"
          visible={visible}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.point_group"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.hall_symbol"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.symmetry.prototype_aflow_id"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuStructure.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuStructure
