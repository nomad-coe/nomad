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
import { screen } from '../../conftest.spec'
import { renderSearchEntry } from '../conftest.spec'
import { WidgetScatterPlotEdit } from './WidgetScatterPlotEdit'

describe('test edit dialog error messages', () => {
  test.each([
    ['missing x', {x: {quantity: 'results.material.n_elements'}}, 'Please specify a value'],
    ['missing y', {y: {quantity: 'results.material.n_elements'}}, 'Please specify a value'],
    ['unavailable x', {x: {quantity: 'results.material.not_a_quantity'}}, 'The quantity "results.material.not_a_quantity" is not available.'],
    ['unavailable y', {y: {quantity: 'results.material.not_a_quantity'}}, 'The quantity "results.material.not_a_quantity" is not available.'],
    ['unavailable color', {markers: {color: {quantity: 'results.material.not_a_quantity'}}}, 'The quantity "results.material.not_a_quantity" is not available.'],
    ['invalid jmespath x', {x: {quantity: 'results.material.n_elements[*'}}, 'Invalid JMESPath query, please check your syntax.'],
    ['invalid jmespath y', {y: {quantity: 'results.material.n_elements[*'}}, 'Invalid JMESPath query, please check your syntax.'],
    ['invalid jmespath color', {markers: {color: {quantity: 'results.material.n_elements[*'}}}, 'Invalid JMESPath query, please check your syntax.'],
    ['no jmespath for repeating x', {x: {quantity: 'results.material.topology.cell.a'}}, 'The quantity "results.material.topology.cell.a" is contained in at least one repeatable section. Please use JMESPath syntax to select one or more target sections.'],
    ['no jmespath for repeating y', {y: {quantity: 'results.material.topology.cell.a'}}, 'The quantity "results.material.topology.cell.a" is contained in at least one repeatable section. Please use JMESPath syntax to select one or more target sections.'],
    ['no jmespath for repeating color', {markers: {color: {quantity: 'results.material.topology.cell.a'}}}, 'The quantity "results.material.topology.cell.a" is contained in at least one repeatable section. Please use JMESPath syntax to select one or more target sections.'],
    ['invalid x unit', {x: {quantity: 'results.material.topology[0].cell.a', unit: 'nounit'}}, 'Unit "nounit" not found.'],
    ['invalid y unit', {y: {quantity: 'results.material.topology[0].cell.a', unit: 'nounit'}}, 'Unit "nounit" not found.'],
    ['invalid color unit', {markers: {color: {quantity: 'results.material.topology[0].cell.a', unit: 'nounit'}}}, 'Unit "nounit" not found.'],
    ['incompatible x unit', {x: {quantity: 'results.material.topology[0].cell.a', unit: 'joule'}}, 'Unit "joule" is incompatible with dimension "length"'],
    ['incompatible y unit', {y: {quantity: 'results.material.topology[0].cell.a', unit: 'joule'}}, 'Unit "joule" is incompatible with dimension "length"'],
    ['incompatible color unit', {markers: {color: {quantity: 'results.material.topology[0].cell.a', unit: 'joule'}}}, 'Unit "joule" is incompatible with dimension "length"']
  ])('%s', async (name, config, error) => {
    const finalConfig = {
      id: '0',
      editing: true,
      ...config
    }
    renderSearchEntry(<WidgetScatterPlotEdit widget={finalConfig} />)
    const button = screen.getByText('Done')
    await userEvent.click(button)
    screen.getByText(error)
  })
})
