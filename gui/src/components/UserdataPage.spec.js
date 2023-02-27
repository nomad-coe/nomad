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
import { render, startAPI, closeAPI } from './conftest.spec'
import { expectFilterMainMenu, expectSearchResults } from './search/conftest.spec'
import { ui } from '../config'
import UserDatapage from './UserdataPage'
import { minutes } from '../setupTests'

test('renders user data search page correctly', async () => {
  const context = ui.apps.options.entries
  await startAPI('tests.states.search.search', 'tests/data/search/userdatapage', 'test', 'password')
  render(<UserDatapage />)

  await expectFilterMainMenu(context)
  await expectSearchResults(context)
  closeAPI()
}, 5 * minutes)
