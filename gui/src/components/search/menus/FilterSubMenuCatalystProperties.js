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
import InputSection from '../input/InputSection'
import InputRange from '../input/InputRange'
import InputField from '../input/InputField'

const FilterSubMenuCatalyst = React.memo(({
  id,
  ...rest
}) => {
  const {selected, open} = useContext(filterMenuContext)
  const visible = open && id === selected

  return <FilterSubMenu id={id} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.properties.catalytic.reaction.name"
          visible={visible}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.catalytic.reaction.reactants"
          visible={visible}
        >
          <InputField
          quantity="results.properties.catalytic.reaction.reactants.name"
          visible={visible}
          />
          <InputRange
          quantity="results.properties.catalytic.reaction.reactants.conversion"
          visible={visible}
          />
          <InputRange
          quantity="results.properties.catalytic.reaction.reactants.gas_concentration_in"
          visible={visible}
          />
          <InputRange
          quantity="results.properties.catalytic.reaction.reactants.gas_concentration_out"
          visible={visible}
          />
        </InputSection>
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.catalytic.reaction.products"
          visible={visible}
        >
          <InputField
          quantity="results.properties.catalytic.reaction.products.name"
          visible={visible}
          />
          <InputRange
          quantity="results.properties.catalytic.reaction.products.selectivity"
          visible={visible}
          />
          <InputRange
          quantity="results.properties.catalytic.reaction.products.gas_concentration_out"
          visible={visible}
          />
        </InputSection>
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputRange
          quantity="results.properties.catalytic.reaction.temperature"
          visible={visible}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.catalytic.catalyst_synthesis"
          disableHeader
          visible={visible}
        >
          <InputField
            quantity="results.properties.catalytic.catalyst_synthesis.catalyst_type"
            visible={visible}
          />
          <InputField
            quantity="results.properties.catalytic.catalyst_synthesis.preparation_method"
            visible={visible}
          />
          <InputField
            quantity="results.properties.catalytic.catalyst_synthesis.catalyst_name"
            visible={visible}
          />
        </InputSection>
        <InputSection
          section="results.properties.catalytic.catalyst_characterization"
          visible={visible}
        >
          <InputField
            quantity="results.properties.catalytic.catalyst_characterization.method"
            visible={visible}
          />
          <InputRange
            quantity="results.properties.catalytic.catalyst_characterization.surface_area"
            visible={visible}
          />
        </InputSection>
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuCatalyst.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuCatalyst
