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
import { startAPI, closeAPI } from '../../conftest.spec'
import { renderSearchEntry, expectInputList } from '../conftest.spec'
import InputList from './InputList'

// Mock the useResizeObserver hook. The test environment does not provide any
// resize events on it's own.
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    // Size is adjusted so that 3 items can fit.
    return {height: 120, ref: undefined}
  }}
})

describe('initial state is loaded correctly', () => {
  beforeAll(async () => {
    await startAPI('tests.states.search.search', 'tests/data/search/inputlist')
  })
  afterAll(() => closeAPI())

  test.each([
    [
      'show all',
      'results.method.simulation.program_name',
      ['VASP', 'exciting'],
      'all'
    ],
    [
      'show n first',
      'results.method.simulation.dft.xc_functional_names',
      ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL', 'LDA_C_PZ'],
      'top'
    ],
    [
      'show only item',
      'results.method.simulation.dft.xc_functional_type',
      ['not processed'],
      'single'
    ],
    [
      'no results',
      'results.material.functional_type',
      [],
      undefined
    ]
  ])('%s', async (name, quantity, items, prompt) => {
    renderSearchEntry(<InputList
      quantity={quantity}
      visible
      draggable
      aggId="statistics"
    />
    )
    await expectInputList(quantity, false, items, prompt)
  })
})
