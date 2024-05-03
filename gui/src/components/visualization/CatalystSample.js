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
 * Shows a summary of catalyst properties.
 */
const CatalystSample = React.memo(({data}) => {
  const prefix = 'results.properties.catalytic.catalyst_characterization'
  const prefix1 = 'results.properties.catalytic.catalyst_synthesis'
  const prefix2 = 'results.properties.catalytic.reaction'
  return data !== false
    ? data
      ? <QuantityTable>
        <QuantityRow>
          <QuantityCell value={data?.catalyst_name} quantity={`${prefix}.catalyst_name`}/>
          <QuantityCell value={data?.surface_area} quantity={`${prefix}.surface_area`}/>
          <QuantityCell value={data?.catalyst_type} quantity={`${prefix1}.catalyst_type`}/>
          <QuantityCell value={data?.preparation_method} quantity={`${prefix1}.preparation_method`}/>
        </QuantityRow>
        <QuantityRow>
          <QuantityCell value={data?.name} quantity={`${prefix2}.name`}/>
          <QuantityCell value={data?.temperature} quantity={`${prefix2}.temperature`}/>
          <QuantityCell value={data?.reactants.name} quantity={`${prefix2}.reactants.name`}/>
          <QuantityCell value={data?.reactants.conversion} quantity={`${prefix2}.reactants.conversion`}/>
          <QuantityCell value={data?.products.name} quantity={`${prefix2}.products.name`}/>
          <QuantityCell value={data?.products.selectivity} quantity={`${prefix2}.products.selectivity`}/>

        </QuantityRow>
      </QuantityTable>
      : <Placeholder />
    : <NoData />
})

CatalystSample.propTypes = {
  data: PropTypes.any
}

export default CatalystSample
