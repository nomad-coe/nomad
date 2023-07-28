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
import { expectPlot, VisualizationState } from './conftest.spec'
import Spectra, { spectraError } from './Spectra'

test.each([
  ['no data', VisualizationState.NoData, false, undefined],
  ['loading', VisualizationState.Loading, undefined, 'spectra-placeholder'],
  ['error: invalid data layout', VisualizationState.Error, {invalid: 'data'}, undefined],
  ['valid', VisualizationState.Success, [{energies: [0, 1, 2], intensities: [0, 2, 1], label: 'computation', type: 'XAS'}], undefined]
])('spectra: %s', async (id, state, data, spectraID) => {
  renderNoAPI(<Spectra data={data} />)
  await expectPlot(state, spectraID, spectraError)
})
