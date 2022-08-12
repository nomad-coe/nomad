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
import assert from 'assert'
import { within } from '@testing-library/dom'
import { render } from '../conftest.spec'
import GeometryOptimization from './GeometryOptimization'
import { PlotState } from './Plot'

  test.each([
    [PlotState.NoData, false],
    [PlotState.Loading, undefined],
    [PlotState.Error, {invalid: "data"}],
    [PlotState.Success, [0, 1, 2, 3, 4]]
  ])('%s', async (state, energies) => {
    const {container} = render(<GeometryOptimization energies={energies}/>)
    await expectGeometryOptimization(state, container)
  })

/**
 * Tests that the geometry optimization data is loaded properly.
 *
 * @param {PlotState} state The expected state of the component.
 * @param {object} root The root element to perform the search on.
 */
async function expectGeometryOptimization(state, container = document.body) {
  assert(state in PlotState, 'Please provide a valid state.')
  const root = within(container)

  if (state === PlotState.NoData) {
    // The component should immediately (without any placeholders) display that
    // there is no data.
    expect(root.getByText('no data')).toBeInTheDocument()
  } else if (state === PlotState.Loading) {
    // The component should immediately display the placeholder
    expect(root.getByTestId('geometry-optimization-plot-placeholder')).toBeInTheDocument()
  } else if (state === PlotState.Error) {
    // The component should immediately (without any placeholders) display the
    // error message.
    expect(root.getByText('Could not load geometry optimization data.')).toBeInTheDocument()
  } else if (state === PlotState.Success) {
    // The component should display the plot. Only way to currently test SVG
    // charts is to use querySelector which is not ideal. Note that testing for
    // SVG text contents in JSDOM is also not currently working with JSDOM. Here
    // we simply test that the SVG elements holding the titles exist.
    expect(container.querySelector(".xtitle")).toBeInTheDocument()
    expect(container.querySelector(".ytitle")).toBeInTheDocument()
    expect(container.querySelector(".y2title")).toBeInTheDocument()
  }
}
