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
import { render, screen } from './conftest.spec'
import {CodeList} from './About'

test('list of codes renders correctly', async () => {
  render(<CodeList />)
  // Check that the list of parsers is printed correctly within the correct categories.
  // This list should be updated if new codes are added or the code names are updated.
  // This test value is hardcoded so that any unintentional changes to the code names and
  // categories can be avoided.
  const list = /Atomistic codes: VASP/
  expect(screen.getByTestId('code-list')).toHaveTextContent(list)
})
