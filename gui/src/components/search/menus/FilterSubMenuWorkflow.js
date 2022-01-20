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
import InputSlider from '../input/InputSlider'
import InputSection from '../input/InputSection'
import InputField from '../input/InputField'
import { useUnits } from '../../../units'

const FilterSubMenuWorkflow = React.memo(({
  value,
  ...rest
}) => {
  const units = useUnits()
  const {selected} = useContext(filterMenuContext)
  const visible = value === selected

  return <FilterSubMenu value={value} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputField
          quantity="workflow"
          visible={visible}
          disableSearch
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.geometry_optimization"
          visible={visible}
        >
          <InputSlider
            quantity="results.properties.geometry_optimization.final_energy_difference"
            units={units}
            visible={visible}
          />
          <InputSlider
            quantity="results.properties.geometry_optimization.final_force_maximum"
            units={units}
            visible={visible}
          />
          <InputSlider
            quantity="results.properties.geometry_optimization.final_displacement_maximum"
            units={units}
            visible={visible}
          />
        </InputSection>
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuWorkflow.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuWorkflow
