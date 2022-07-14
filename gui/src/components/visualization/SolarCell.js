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
import React, { } from 'react'
import PropTypes from 'prop-types'
import { QuantityTable, QuantityRow, QuantityCell } from '../Quantity'
import NoData from './NoData'
import Placeholder from './Placeholder'

/**
 * Shows a summary of solar cell properties.
 */
const SolarCell = React.memo(({data}) => {
  const prefix = 'results.properties.optoelectronic.solar_cell'
  return data !== false
    ? data
      ? <QuantityTable>
        <QuantityRow>
          <QuantityCell value={data?.efficiency} quantity={`${prefix}.efficiency`}/>
          <QuantityCell value={data?.fill_factor} quantity={`${prefix}.fill_factor`}/>
          <QuantityCell value={data?.open_circuit_voltage} quantity={`${prefix}.open_circuit_voltage`}/>
          <QuantityCell value={data?.short_circuit_current_density} quantity={`${prefix}.short_circuit_current_density`}/>
        </QuantityRow>
        <QuantityRow>
          <QuantityCell value={data?.illumination_intensity} quantity={`${prefix}.illumination_intensity`}/>
          <QuantityCell value={data?.device_area} quantity={`${prefix}.device_area`}/>
          <QuantityCell value={data?.device_architecture} colSpan={2} quantity={`${prefix}.device_architecture`}/>
        </QuantityRow>
        <QuantityRow>
          <QuantityCell value={data?.absorber} colSpan={2} quantity={`${prefix}.absorber`}/>
          <QuantityCell value={data?.absorber_fabrication} colSpan={2} quantity={`${prefix}.absorber_fabrication`}/>
        </QuantityRow>
        <QuantityRow>
          <QuantityCell value={data?.device_stack} colSpan={4} maxWidth="none" quantity={`${prefix}.device_stack`}/>
        </QuantityRow>
      </QuantityTable>
      : <Placeholder />
    : <NoData />
})

SolarCell.propTypes = {
  data: PropTypes.any
}

export default SolarCell
