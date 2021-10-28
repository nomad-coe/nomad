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
import InputSelect from '../input/InputSelect'
import InputSlider from '../input/InputSlider'
import InputCheckboxes from '../input/InputCheckboxes'
import InputSection from '../input/InputSection'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import { Quantity, useUnits } from '../../../units'

const step = new Quantity(10, 'gigapascal')
const FilterSubMenuElectronic = React.memo(({
  value,
  ...rest
}) => {
  const {selected} = useContext(filterMenuContext)
  const visible = value === selected
  const units = useUnits()

  return <FilterSubMenu value={value} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.mechanical.bulk_modulus"
          visible={visible}
        >
          <InputSelect
            quantity="results.properties.mechanical.bulk_modulus.type"
            visible={visible}
          />
          <InputSlider
            quantity="results.properties.mechanical.bulk_modulus.value"
            units={units}
            step={step}
            visible={visible}
          />
        </InputSection>
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.mechanical.shear_modulus"
          visible={visible}
        >
          <InputSelect
            quantity="results.properties.mechanical.shear_modulus.type"
            visible={visible}
          />
          <InputSlider
            quantity="results.properties.mechanical.shear_modulus.value"
            units={units}
            step={step}
            visible={visible}
          />
        </InputSection>
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputCheckboxes
          quantity="mechanical_properties"
          visible={visible}
        />
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuElectronic.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuElectronic
