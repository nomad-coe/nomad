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
import { WidgetHistogramEdit } from './WidgetHistogramEdit'

describe('test edit dialog error messages', () => {
  test.each([
    ['missing x', {x: {}}, 'Please specify a value.'],
    ['unavailable x', {x: {search_quantity: 'results.material.not_a_quantity'}}, 'The quantity "results.material.not_a_quantity" is not available.'],
    ['invalid x unit', {x: {search_quantity: 'results.material.topology.cell.a', unit: 'nounit'}}, 'Unit "nounit" not found.'],
    ['incompatible x unit', {x: {search_quantity: 'results.material.topology.cell.a', unit: 'joule'}}, 'Unit "joule" is incompatible with dimension "length".']
  ])('%s', async (name, config, error) => {
    const finalConfig = {
      id: '0',
      editing: true,
      ...config
    }
    renderSearchEntry(<WidgetHistogramEdit widget={finalConfig} />)
    const button = screen.getByText('Done')
    await userEvent.click(button)
    screen.getByText(error)
  })
})
