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
import InputRange from '../input/InputRange'
import InputField from '../input/InputField'
import InputSection from '../input/InputSection'
import { InputGrid, InputGridItem } from '../input/InputGrid'

const FilterSubMenuMolecularDynamics = React.memo(({
  id,
  ...rest
}) => {
  const {selected} = useContext(filterMenuContext)
  const visible = id === selected

  return <FilterSubMenu id={id} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.thermodynamic.trajectory"
          disableHeader
          visible={visible}
        >
          <InputField
            quantity="results.properties.thermodynamic.trajectory.available_properties"
            visible={visible}
            disableSearch
            formatLabels
          />
          <InputField
            quantity="results.properties.thermodynamic.trajectory.provenance.molecular_dynamics.ensemble_type"
            disableSearch
            visible={visible}
          />
          <InputRange
            quantity="results.properties.thermodynamic.trajectory.provenance.molecular_dynamics.time_step"
            visible={visible}
          />
        </InputSection>
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuMolecularDynamics.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuMolecularDynamics
