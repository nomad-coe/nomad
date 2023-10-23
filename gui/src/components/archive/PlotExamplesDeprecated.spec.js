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
import {render, screen, within} from '../conftest.spec'
import {PlotExamplesDeprecated} from './PlotExamplesDeprecated'

test('correctly renders the deprecated plots', async () => {
  render(<PlotExamplesDeprecated />)

  await screen.findByText('Examples to show how to plot simple and multiline graphs and how to customize the line styles')
  await screen.findByText('Examples to plot data from an array of section')

  const plots = screen.queryAllByText('plot')
  expect(plots.length).toBe(8)
  const plotCards = plots.map(plot => plot.parentElement.parentElement) // get the card

  within(plotCards[0]).getByText('Process Time (fs)')
  within(plotCards[0]).getByText('T Substrate (K)')

  within(plotCards[1]).getByText('Process Time (fs)')
  within(plotCards[1]).getByText('Temperature (K)')
  within(plotCards[1]).getByText('Set Substrate Temperature')
  within(plotCards[1]).getByText('Substrate Temperature')

  within(plotCards[2]).getByText('Process Time (fs)')
  within(plotCards[2]).getByText('Substrate Temperature (K)')
  within(plotCards[2]).getByText('Chamber Pressure (GPa)')

  within(plotCards[3]).getByText('Process Time (fs)')
  within(plotCards[3]).getByText('Substrate Temperature (K)')
  within(plotCards[3]).getByText('Chamber Pressure (GPa)')

  within(plotCards[4]).getByText('Process Time (fs)')
  within(plotCards[4]).getByText('Chamber Pressure (GPa)')

  within(plotCards[5]).getByText('Process Time (fs)')
  within(plotCards[5]).getByText('Substrate Temperature (K)')

  // Check that first subsection is plotted even if third does not exist
  within(plotCards[6]).getByText('Process Time (fs)')
  within(plotCards[6]).getByText('Substrate Temperature (K)')

  // Check that legend for repeating subsection gets subsection name
  within(plotCards[7]).getByText('Substrate 1, Substrate Temperature')
  within(plotCards[7]).getByText('Substrate 2, Substrate Temperature')
})
