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
import { ui } from '../../config'
import { render, screen } from '../conftest.spec'
import FilterSummary from './FilterSummary'
import { SearchContext } from './SearchContext'

test.each([
  ['integer', 'results.material.n_elements', 12, 'N elements', 'n_elements=12'],
  ['string', 'results.material.symmetry.crystal_system', 'cubic', 'Crystal system', 'cubic'],
  ['float', 'results.method.simulation.precision.k_line_density', 12.3, 'k-line density (Å)', 'k_line_density=12.3'],
  ['datetime', 'upload_create_time', 0, 'Upload create time', 'upload_create_time=01/01/1970'],
  ['boolean', 'results.properties.electronic.dos_electronic.spin_polarized', 'false', 'Spin-polarized', 'false']
])('%s', async (name, quantity, input, title, output) => {
  const context = ui.apps.options.entries
  render(
    <SearchContext
        resource={context.resource}
        initialPagination={context.pagination}
        initialColumns={context.columns}
        initialRows={context.rows}
        initialFilterMenus={context.filter_menus}
        initialFiltersLocked={context.filters_locked}
        initialDashboard={context?.dashboard}
        initialFilterValues={{[quantity]: input}}
        initialSearchSyntaxes={context?.search_syntaxes}
    >
      <FilterSummary quantities={new Set([quantity])}/>
    </SearchContext>
  )
  expect(screen.getByText(title, {exact: false})).toBeInTheDocument()
  expect(screen.getByText(output)).toBeInTheDocument()
})
