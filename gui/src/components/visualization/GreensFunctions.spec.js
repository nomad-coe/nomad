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
import { expectPlot, VisualizationState } from './conftest.spec'
import GreensFunctions, { gfError } from './GreensFunctions'

test.each([
  ['no data', VisualizationState.NoData, false, undefined, undefined, undefined, undefined],
  ['loading', VisualizationState.Loading, undefined, undefined, undefined, 'greens-functions-regtau', 'greens-functions-imsiw'],
  ['error: invalid data layout', VisualizationState.Error, {invalid: 'data'}, undefined, undefined, undefined, undefined],
  ['valid', VisualizationState.Success, {tau: [0, 1, 2], regtau: [[[[1, 0, 1]], [[1, 0, 1]]]], iw: [-1, 0, 1], imsiw: [[[[1, 0, 1]], [[1, 0, 1]]]]}, {magnetic_state: 'paramagnetic'}, undefined, 'greens-functions-regtau', 'greens-functions-imsiw']
])('greens functions: %s', async (id, state, data, methodology, classes, regtauTestID, imsiwTestID) => {
  render(<GreensFunctions data={data} methodology={methodology} classes={classes} />)
  await expectPlot(state, regtauTestID, gfError)
  await expectPlot(state, imsiwTestID, gfError)
})
