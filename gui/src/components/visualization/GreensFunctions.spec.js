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
import GreensFunctions, { gfError } from './GreensFunctions'

const greens_function_freq_im = 'gfs-greens_function_freq_im-placeholder'

test.each([
  ['no data', VisualizationState.NoData, {greens_function_freq_im: false}, undefined],
  ['loading', VisualizationState.Loading, {greens_function_freq_im: undefined}, undefined],
  ['error: invalid data layout', VisualizationState.Error, {invalid: 'data'}, undefined],
  [
    'valid',
    VisualizationState.Success,
    {greens_function_freq_im: [{x: [0, 1], y: [[0, 1]]}]},
    {magnetic_state: 'paramagnetic'}
  ]
])('greens functions: %s', async (id, state, data, provenance) => {
  renderNoAPI(<GreensFunctions data={data} provenance={provenance} />)
  await expectPlot(state, greens_function_freq_im, gfError)
})
