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
import React from 'react'
import { render } from '../conftest.spec'
import { expectMethodologyItem } from '../entry/conftest.spec'
import { expectPlot } from './conftest.spec'
import { PlotState } from './Plot'
import RadialDistributionFunction, { rdfError, rdfPath } from './RadialDistributionFunction'

test.each([
  ['no data', PlotState.NoData, {molecular: {'MOL-MOL': false}}, undefined],
  ['loading', PlotState.Loading, {molecular: {'MOL-MOL': undefined}}, 'radial-distribution-function-molecular-mol-mol-placeholder'],
  ['error: data cannot be false', PlotState.Error, false, undefined],
  ['error: data cannot be undefined', PlotState.Error, undefined, undefined],
  ['error: invalid data layout', PlotState.Error, {invalid: "data"}, undefined],
  ['valid', PlotState.Success, {molecular: {'MOL-MOL': [{bins: [0, 1], value: [0, 1]}]}}, undefined]
])('rdf plot: %s', async (id, state, data, placeholderTestID) => {
  render(<RadialDistributionFunction rdf={data} />)
  await expectPlot(state, placeholderTestID, rdfError)
})

test.each([
  ['no methodology', undefined],
  ['valid methodology', {molecular_dynamics: {time_step: 2e-15, ensemble_type: 'NVT'}}]
])('methodology is displayed correctly: %s', async (state, methodology) => {
  render(<RadialDistributionFunction rdf={{}} methodology={methodology}/>)
  expectMethodologyItem(
    'Molecular dynamics',
    methodology,
    `${rdfPath.join('.')}.methodology`
  )
})
