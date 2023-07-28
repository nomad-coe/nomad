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
import InputSection from '../input/InputSection'
import { InputCheckboxValue } from '../input/InputCheckbox'

const FilterSubMenuEELS = React.memo(({
  id,
  ...rest
}) => {
  const {selected, open} = useContext(filterMenuContext)
  const visible = open && id === selected

  return <FilterSubMenu
    id={id}
    actions={<InputCheckboxValue
      quantity="results.method.method_name"
      value="EELS"
      description="Search EELS entries"
    />}
    {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputSection
          section="results.properties.spectroscopic.spectra.provenance.eels"
          disableHeader
          visible={visible}
        >
          <InputField
            quantity="results.properties.spectroscopic.spectra.provenance.eels.detector_type"
            visible={visible}
            xs={12}
          />
          <InputRange
            quantity="results.properties.spectroscopic.spectra.provenance.eels.resolution"
            visible={visible}
          />
          <InputRange
            quantity="results.properties.spectroscopic.spectra.provenance.eels.min_energy"
            visible={visible}
          />
          <InputRange
            quantity="results.properties.spectroscopic.spectra.provenance.eels.max_energy"
            visible={visible}
          />
        </InputSection>
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuEELS.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuEELS
