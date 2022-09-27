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
import GeometryOptimization from './GeometryOptimization'

const errorMsg = 'Could not load geometry optimization data.'
test.each([
  [VisualizationState.NoData, false, undefined],
  [VisualizationState.Loading, undefined, 'geometry-optimization-plot-placeholder'],
  [VisualizationState.Error, {invalid: 'data'}, undefined],
  [VisualizationState.Success, [0, 1, 2, 3, 4], undefined]
])('%s', async (state, energies, placeholderTestID) => {
  const {container} = render(<GeometryOptimization energies={energies}/>)
  await expectPlot(state, placeholderTestID, errorMsg, container)
})
