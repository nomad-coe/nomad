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
import { render, screen } from '../../conftest.spec'
import { InputJMESPath } from './InputMetainfo'
import { SearchContext } from '../SearchContext'
import { ui } from '../../../config'
import { DType } from '../../../utils'

const context = ui.apps.options.entries

test.each([
  ['filter exists', 'results.material.n_elements', null],
  ['filter does not exist', 'not.present', "The quantity \"not.present\" is not available."],
  ['valid jmespath', 'results.material.n_elements[0]', null],
  ['invalid jmespath', 'results.material.n_elements[0', 'Invalid JMESPath query, please check your syntax.']
])('%s', async (name, input, error) => {
  const onErrorMock = jest.fn()
  const onAcceptMock = jest.fn()
  const onChangeMock = jest.fn()
  render(
    <SearchContext
        resource={context.resource}
        initialPagination={context.pagination}
        initialColumns={context.columns}
        initialRows={context.rows}
        initialFilters={context?.filters}
        initialFilterMenus={context.filter_menus}
        initialFiltersLocked={context.filters_locked}
        initialDashboard={context?.dashboard}
        initialSearchSyntaxes={context?.search_syntaxes}
    >
        <InputJMESPath
          value={input}
          onAccept={onAcceptMock}
          onError={onErrorMock}
          onChange={onChangeMock}
          dtypes={new Set([DType.String, DType.Int])}
        />
    </SearchContext>
  )
  const textInput = screen.getByRole('textbox')
  await userEvent.type(textInput, '{Enter}')
  if (error) {
    expect(onErrorMock).toHaveBeenCalledWith(error)
  } else {
    expect(onAcceptMock).toHaveBeenCalledWith(input, undefined)
  }
})
