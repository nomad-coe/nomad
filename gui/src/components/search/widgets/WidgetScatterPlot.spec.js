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
import React, { useMemo } from 'react'
import { screen } from '../../conftest.spec'
import { expectWidgetScatterPlot, renderSearchEntry } from '../conftest.spec'
import { WidgetScatterPlot } from './WidgetScatterPlot'

// Mock the resize observer that defines the widget size
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {height: 400, width: 400, ref: undefined}
  }}
})

// Mock the useHits hook that returns the plotted data. Note that the useHits
// hook needs to memo things to prevent rendering loops.
const mockUseMemo = useMemo
jest.mock('../SearchContext', () => ({
    ...jest.requireActual('../SearchContext'),
    useSearchContext: () => ({
      ...jest.requireActual('../SearchContext').useSearchContext(),
      useHits: (id, required, pagination, callback) => {
        const response = mockUseMemo(() => {
          callback()
          return [
            {
              entry_id: 0,
              entry_create_time: '2024-04-18T09:00:00.000000+00:00',
              results: {
                material: {n_elements: 1, elements: ['Si'], chemical_formula_hill: 'Si2'},
                properties: {
                  electronic: {band_gap: [{value: 0}, {value: 1}]},
                  catalytic: {
                    reaction: {
                      temperature: [0, 1, 2],
                      reactants: [
                        {name: 'CO2', conversion: [0, 1, 2]},
                        {name: 'NO2', conversion: [3, 4, 5]}
                      ]
                    }
                  }
                }
              }
            },
            {
              entry_id: 1,
              entry_create_time: '2024-04-18T10:00:00.000000+00:00',
              results: {
                material: {n_elements: 1, elements: ['C', 'O'], chemical_formula_hill: 'CO2'},
                properties: {
                  electronic: {band_gap: [{value: 0}, {value: 1}]},
                  catalytic: {
                    reaction: {
                      temperature: [0, 1, 2],
                      reactants: [
                        {name: 'H2O', conversion: [0, 1, 2]},
                        {name: 'O2', conversion: [3, 4, 5]}
                      ]
                    }
                  }
                }
              }
            }
          ]
        }, [])
        return response
      }
    })
}))

describe('test different combinations of x/y/color produced with JMESPath', () => {
  test.each([
    ['scalar quantity', 'results.material.n_elements', 'results.material.n_elements', 'results.material.n_elements'],
    ['datetime', 'entry_create_time', 'entry_create_time', 'entry_id'],
    ['index expression', 'results.properties.electronic.band_gap[0].value', 'results.properties.electronic.band_gap[0].value', 'results.properties.electronic.band_gap[0].value'],
    ['slicing', 'results.properties.electronic.band_gap[1:2].value', 'results.properties.electronic.band_gap[1:2].value', 'results.properties.electronic.band_gap[1:2].value'],
    ['function', 'min(results.properties.electronic.band_gap[*].value)', 'min(results.properties.electronic.band_gap[*].value)', 'min(results.properties.electronic.band_gap[*].value)'],
    ['filter projection', "results.material.topology[?label=='original'].cell.a", "results.material.topology[?label=='original'].cell.a", "results.material.topology[?label=='original'].cell.a"],
    ['1D array vs 2D array vs 1D color', 'results.properties.catalytic.reaction.temperature', 'results.properties.catalytic.reaction.reactants[*].conversion', 'results.properties.catalytic.reaction.reactants[*].name']
  ])('%s', async (name, x, y, color) => {
    const config = {
      id: '0',
      scale: 'linear',
      x: {quantity: x},
      y: {quantity: y},
      markers: {color: {quantity: color}}
    }
    renderSearchEntry(<WidgetScatterPlot {...config} />)
    await expectWidgetScatterPlot(config, false)
  })
})

describe('test custom axis titles', () => {
  test.each([
    ['x, no unit', {x: {title: 'My Title', quantity: 'results.material.n_elements'}}, 'My Title'],
    ['y, no unit', {y: {title: 'My Title', quantity: 'results.material.n_elements'}}, 'My Title'],
    ['color, no unit', {markers: {color: {title: 'My Title', quantity: 'results.material.n_elements'}}}, 'My Title'],
    ['x, with unit', {x: {title: 'My Title', quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'My Title (eV)'],
    ['y, with unit', {y: {title: 'My Title', quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'My Title (eV)'],
    ['color, with unit', {markers: {color: {title: 'My Title', quantity: 'results.properties.geometry_optimization.final_energy_difference'}}}, 'My Title (eV)']
  ])('%s', async (name, config, title) => {
    const configFinal = {
      id: '0',
      scale: 'linear',
      x: {quantity: 'results.material.n_elements'},
      y: {quantity: 'results.material.n_elements'},
      markers: {color: {quantity: 'results.material.n_elements'}},
      ...config
    }
    renderSearchEntry(<WidgetScatterPlot {...configFinal} />)
    screen.getByText(title)
  })
})

describe('test custom axis units', () => {
  test.each([
    ['x', {x: {unit: 'Ha', quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'Final energy difference (Ha)'],
    ['y', {y: {unit: 'Ha', quantity: 'results.properties.geometry_optimization.final_energy_difference'}}, 'Final energy difference (Ha)'],
    ['color', {markers: {color: {unit: 'Ha', quantity: 'results.properties.geometry_optimization.final_energy_difference'}}}, 'Final energy difference (Ha)']
  ])('%s', async (name, config, title) => {
    const configFinal = {
      id: '0',
      scale: 'linear',
      x: {quantity: 'results.material.n_elements'},
      y: {quantity: 'results.material.n_elements'},
      markers: {color: {quantity: 'results.material.n_elements'}},
      ...config
    }
    renderSearchEntry(<WidgetScatterPlot {...configFinal} />)
    screen.getByText(title)
  })
})

describe('test different colors', () => {
  test.each([
    ['empty', undefined, undefined, []],
    ['scalar integer', 'results.material.n_elements', {quantity: 'results.material.n_elements'}, []],
    ['scalar float', 'results.properties.geometry_optimization.final_energy_difference', {quantity: 'results.properties.geometry_optimization.final_energy_difference'}, []],
    ['scalar string', 'results.material.chemical_formula_hill', undefined, ['Si2', 'CO2']],
    ['array string', 'results.material.elements', undefined, ['Si', 'C, O']]
  ])('%s', async (name, axis, colorTitle, legend) => {
    const config = {
      id: '0',
      scale: 'linear',
      x: {quantity: 'results.material.n_elements'},
      y: {quantity: 'results.material.n_elements'},
      markers: {color: {quantity: axis}}
    }
    renderSearchEntry(<WidgetScatterPlot {...config} />)
    await expectWidgetScatterPlot(config, false, colorTitle, legend)
  })
})

describe('test error messages', () => {
  test.each([
    ['invalid JMESPath', 'results.properties.electronic.band_gap[*.value', 'Invalid JMESPath query, please check your syntax.']
  ])('%s', async (name, axis, message) => {
    const config = {
      id: '0',
      scale: 'linear',
      x: {quantity: axis},
      y: {quantity: axis},
      markers: {color: {quantity: axis}}
    }
    renderSearchEntry(<WidgetScatterPlot {...config} />)
    screen.getByText(message, {exact: false})
  })
})
