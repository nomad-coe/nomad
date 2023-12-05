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
import { render, screen } from '../conftest.spec'
import { withFilters } from './FilterRegistry'

const TestComponent = withFilters(({initialFilterData}) => {
  return Object.keys(initialFilterData).map((filter) => <div key={filter}>{filter}</div>)
})

test.each([
  ['nomad default', {}, ['mainfile'], []],
  ['nomad include', {include: ['mainfile']}, ['mainfile'], []],
  ['nomad exclude', {exclude: ['mainfile']}, [], ['mainfile']]
])('%s', async (name, config, include, exclude) => {
  render(<TestComponent initialFilters={config}/>)
  for (const filter of include) {
    expect(await screen.findByText(filter, {exact: true})).toBeInTheDocument()
  }
  for (const filter of exclude) {
    await expect(screen.findByText(filter, {exact: true})).rejects.toThrow()
  }
})
