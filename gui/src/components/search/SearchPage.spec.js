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
import { render, startAPI, closeAPI } from '../conftest.spec'
import { expectFilterMainMenu, expectSearchResults } from './conftest.spec'
import { ui } from '../../config'
import { SearchContext } from './SearchContext'
import SearchPage from './SearchPage'
import { minutes } from '../../setupTests'

describe('', () => {
  beforeAll(async () => {
    await startAPI('tests.states.search.search', 'tests/data/search/searchpage')
  })
  afterAll(() => closeAPI())

  test.each(
    Object.entries(ui.apps.options)
  )('renders search page correctly, context: %s', async (key, context) => {
    render(
      <SearchContext
          resource={context.resource}
          initialPagination={context.pagination}
          initialColumns={context.columns}
          initialRows={context.rows}
          initialFilterMenus={context.filter_menus}
          initialFiltersLocked={context.filters_locked}
          initialDashboard={context?.dashboard}
          initialSearchSyntaxes={context?.search_syntaxes}
          id={context?.path}
      >
        <SearchPage />
      </SearchContext>
    )

    await expectFilterMainMenu(context)
    await expectSearchResults(context)
  }, 5 * minutes)
})
