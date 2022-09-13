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

import assert from 'assert'
import { within } from '@testing-library/dom'
import { PlotState } from './Plot'

/*****************************************************************************/
// Expects

/**
 * Tests that all plot buttons are in place.
 *
 * @param {object} root The container to work on.
 */
export function expectPlotButtons(container) {
  const root = within(container)
  expect(root.getByRole('button', {name: 'Reset view'})).toBeInTheDocument()
  expect(root.getByRole('button', {name: 'Toggle fullscreen'})).toBeInTheDocument()
  expect(root.getByRole('button', {name: 'Capture image'})).toBeInTheDocument()
  expect(root.getByRole('button', {name: 'View data in the archive'})).toBeInTheDocument()
}

/**
 * Tests that plot is loaded properly.
 *
 * @param {PlotState} state The expected plot state.
 * @param {str} placeholderTestID The test id for the placeholder.
 * @param {str} errorMsg The expected error message.
 * @param {object} container The root element to perform the search on.
 */
export async function expectPlot(
  state,
  placeholderTestID,
  errorMsg,
  container = document.body
) {
  assert(state in PlotState, 'Please provide a valid state.')
  const root = within(container)

  if (state === PlotState.NoData) {
    // The component should immediately (without any placeholders) display that
    // there is no data.
    expect(root.getByText('no data')).toBeInTheDocument()
  } else if (state === PlotState.Loading) {
    // The component should immediately display the placeholder
    expect(root.getByTestId(placeholderTestID)).toBeInTheDocument()
  } else if (state === PlotState.Error) {
    // The component should immediately (without any placeholders) display the
    // error message.
    expect(root.getByText(errorMsg)).toBeInTheDocument()
  } else if (state === PlotState.Success) {
    // The component should display the plot. Only way to currently test SVG
    // charts is to use querySelector which is not ideal. Note that testing for
    // SVG text contents in JSDOM is also not currently working with JSDOM. Here
    // we simply test that the SVG elements holding the titles exist.
    expect(container.querySelector(".xtitle")).toBeInTheDocument()
    expect(container.querySelector(".ytitle")).toBeInTheDocument()
  }
}
