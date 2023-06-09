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
import InputRange from '../input/InputRange'
import { InputCheckboxValue } from '../input/InputCheckbox'

const FilterSubMenuDFT = React.memo(({
  id,
  ...rest
}) => {
  const {selected, open} = useContext(filterMenuContext)
  const visible = open && id === selected

  return <FilterSubMenu
    id={id}
    actions={<InputCheckboxValue
      quantity="results.method.method_name"
      value="DFT"
      description="Search DFT entries"
    />}
    {...rest}
  >
    <InputGrid>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.method.simulation.dft.xc_functional_type"
          visible={visible}
          xs={12}
          disableSearch
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.method.simulation.dft.xc_functional_names"
          visible={visible}
          xs={12}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputRange
          quantity="results.method.simulation.dft.exact_exchange_mixing_factor"
          visible={visible}
          xs={12}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputRange
          quantity="results.method.simulation.dft.hubbard_kanamori_model.u_effective"
          visible={visible}
          xs={12}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.method.simulation.dft.core_electron_treatment"
          visible={visible}
          disableSearch
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.method.simulation.dft.relativity_method"
          visible={visible}
          disableSearch
        />
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuDFT.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuDFT
