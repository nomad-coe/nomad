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
import userEvent from '@testing-library/user-event'
import { screen, renderNoAPI } from '../conftest.spec'
import DOS from './DOS'

test.each([
  [
    'no normalization factors',
    [{
        energies: [0.1, 0.2, 0.3, 0.4, 0.5],
        densities: [[0.1, 0.2, 0.3, 0.4, 0.5]],
        normalization_factors: undefined
    }],
    true
  ],
  [
    'with normalization factors',
    [{
        energies: [0.1, 0.2, 0.3, 0.4, 0.5],
        densities: [[0.1, 0.2, 0.3, 0.4, 0.5]],
        normalization_factors: [2.5]
    }],
    false
  ]
])('DOS: %s', async (id, data, normalization_disabled) => {
  // Render component
  renderNoAPI(<DOS data={data} />)

  // Find menu element
  const menuButton = screen.getByTestId('dos-options-menu')

  // Click on menu to open options
  await userEvent.click(menuButton)

  // Assert that the normalization option is disabled or not
  const normalizationOption = screen.getByLabelText('Normalize intensities')
  if (normalization_disabled) {
    expect(normalizationOption).toBeDisabled()
  } else {
    expect(normalizationOption).not.toBeDisabled()
  }
})
