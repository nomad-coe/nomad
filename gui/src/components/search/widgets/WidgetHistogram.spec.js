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
import { renderSearchEntry, expectWidgetHistogram } from '../conftest.spec'
import { WidgetHistogram } from './WidgetHistogram'

// Mock the resize observer that defines the widget size
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {height: 400, width: 400, ref: undefined}
  }}
})

describe.only('test custom axis titles', () => {
  test.each([
    ['no unit', {x: {title: 'My Title', search_quantity: 'results.material.n_elements'}}, 'My Title'],
    ['with unit', {x: {title: 'My Title', search_quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'My Title (eV)']
  ])('%s', async (name, config, title) => {
    const configFinal = {
      id: '0',
      y: {scale: 'linear'},
      ...config
    }
    renderSearchEntry(<WidgetHistogram {...configFinal} />)
    await expectWidgetHistogram(configFinal)
  })
})

describe('test custom axis units', () => {
  test.each([
    ['x', {x: {unit: 'Ha', search_quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'Final energy difference (Ha)'],
    ['y', {y: {unit: 'Ha', search_quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'Final energy difference (Ha)'],
    ['color', {markers: {color: {unit: 'Ha', search_quantity: 'results.properties.geometry_optimization.final_energy_difference'}}}, 'Final energy difference (Ha)']
  ])('%s', async (name, config, title) => {
    const configFinal = {
      id: '0',
      scale: 'linear',
      x: {search_quantity: 'results.material.n_elements'},
      y: {search_quantity: 'results.material.n_elements'},
      markers: {color: {search_quantity: 'results.material.n_elements'}},
      ...config
    }
    renderSearchEntry(<WidgetHistogram {...configFinal} />)
    await expectWidgetHistogram(configFinal)
  })
})
