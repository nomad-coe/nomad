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
import InputSlider from '../input/InputSlider'
import { Quantity, useUnits } from '../../../units'

const stepResolution = new Quantity(0.1, 'electron_volt')
const stepEnergyWindow = new Quantity(100, 'electron_volt')

const FilterSubMenuEELS = React.memo(({
  value,
  ...rest
}) => {
  const units = useUnits()
  const {selected} = useContext(filterMenuContext)
  const visible = value === selected

  return <FilterSubMenu value={value} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputSlider
          quantity="results.method.experiment.eels.resolution"
          units={units}
          step={stepResolution}
          visible={visible}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputSlider
          quantity="results.method.experiment.eels.energy_window"
          units={units}
          step={stepEnergyWindow}
          visible={visible}
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.method.experiment.eels.detector_type"
          visible={visible}
          xs={12}
        />
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuEELS.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuEELS
