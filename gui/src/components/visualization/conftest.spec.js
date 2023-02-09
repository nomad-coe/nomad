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
import { within } from '../conftest.spec'
import { webGlError } from '../ErrorHandler'

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
 * Tests that a visualization state is loaded correctly.
 *
 * @param {VisualizationState} state The expected visualization state.
 * @param {str} placeholderTestID The test id for the placeholder.
 * @param {str} errorMsg The expected error message.
 * @param {object} container The root element to perform the search on.
 */
export async function expectVisualization(
  state,
  placeholderTestID,
  errorMsg,
  container = document.body
) {
  assert(state in VisualizationState, 'Please provide a valid state.')
  const root = within(container)

  if (state === VisualizationState.NoData) {
    // The component should immediately (without any placeholders) display that
    // there is no data.
    expect(root.queryAllByText('no data').length).toBeGreaterThan(0)
  } else if (state === VisualizationState.Loading) {
    // The component should immediately display the placeholder
    expect(root.getByTestId(placeholderTestID)).toBeInTheDocument()
  } else if (state === VisualizationState.Error) {
    // The component should immediately (without any placeholders) display the
    // error message.
    expect(root.getByText(errorMsg)).toBeInTheDocument()
  } else if (state === VisualizationState.NoWebGL) {
    // The component should immediately (without any placeholders) display the
    // error message.
    expect(root.getByText(webGlError)).toBeInTheDocument()
  }
}

/**
 * Tests that plot is loaded properly.
 *
 * @param {VisualizationState} state The expected plot state.
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
  assert(state in VisualizationState, 'Please provide a valid state.')
  if (state === VisualizationState.Success) {
    // The component should display the plot. Only way to currently test SVG
    // charts is to use querySelector which is not ideal. Note that testing for
    // SVG text contents in JSDOM is also not currently working with JSDOM. Here
    // we simply test that the SVG elements holding the titles exist.
    expect(container.querySelector(".xtitle")).toBeInTheDocument()
    expect(container.querySelector(".ytitle")).toBeInTheDocument()
  } else {
    await expectVisualization(state, placeholderTestID, errorMsg, container)
  }
}

export const VisualizationState = {
  NoData: 'NoData',
  Loading: 'Loading',
  Success: 'Success',
  Error: 'Error',
  NoWebGL: 'NoWebGL'
}

/**
 * Tests that the Pagination component is shown properly.
 *
 * @param {bool} showMore Is show more button visible
 * @param {bool} showLess Is show less button visible
 * @param {bool} loadingMore Is loading more indicator and tooltip present
 */
export async function expectPagination(
  showMore,
  showLess,
  loadingMore,
  container = document.body
) {
  const root = within(container)
  const moreButton = root.queryByButtonText('Show more')
  const lessButton = root.queryByButtonText('Show less')
  showMore
    ? expect(moreButton).toBeInTheDocument()
    : expect(moreButton).toBe(null)
  showLess
    ? expect(lessButton).toBeInTheDocument()
    : expect(lessButton).toBe(null)
  loadingMore && expect(moreButton).toHaveAttribute('disabled')
}
