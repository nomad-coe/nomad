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
import 'regenerator-runtime/runtime'
import { within } from '@testing-library/react'
import { render, screen } from '../conftest.spec'
import MainMenu from './MainMenu'
import { routes } from './Routes'

describe('<MainMenu />', () => {
  it('renders menu items for all nav paths', async () => {
    render(<MainMenu/>)
    for (const key of Object.keys(routes)) {
      const route = routes[key]
      if (route.menu) {
        const menu = await screen.findByTestId(route.menu)
        expect(menu).toBeInTheDocument()
        if (route.routes) {
          for (const subRoute of route.routes) {
            if (subRoute.menu) {
              const item = await within(menu).findByTestId(subRoute.menu)
              expect(item).toBeInTheDocument()
            }
          }
        }
      }
    }
  })
})
