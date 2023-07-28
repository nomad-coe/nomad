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
import { renderNoAPI } from '../conftest.spec'
import { expectMethodologyItem } from '../entry/conftest.spec'
import { expectPlot, VisualizationState } from './conftest.spec'
import Trajectory, { trajectoryError, trajectoryPath } from './Trajectory'

test.each([
  ['no data', VisualizationState.NoData, false, false, false, undefined],
  ['loading', VisualizationState.Loading, undefined, false, false, 'trajectory-placeholder'],
  ['error', VisualizationState.Error, {invalid: 'data'}, false, false, undefined],
  ['valid', VisualizationState.Success, {time: [0, 1, 2], value: [0, 1, 2]}, false, false, undefined]
])('trajectory plot: %s', async (id, state, temperature, pressure, energyPotential, placeholderTestID) => {
  renderNoAPI(<Trajectory
    temperature={temperature}
    pressure={pressure}
    energyPotential={energyPotential}
  />)
  await expectPlot(state, placeholderTestID, trajectoryError)
})

test.each([
  ['no provenance', undefined],
  ['valid provenance', {molecular_dynamics: {time_step: 2e-15, ensemble_type: 'NVT'}}]
])('provenance is displayed correctly: %s', async (state, provenance) => {
  renderNoAPI(<Trajectory provenance={provenance}/>)
  expectMethodologyItem(
    'Molecular dynamics',
    provenance,
    `${trajectoryPath.join('.')}.provenance`
  )
})
