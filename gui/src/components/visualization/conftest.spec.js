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

import { within } from '@testing-library/react'

/*****************************************************************************/
// Expects

/**
 * Tests that all plot buttons are in place.
 *
 * @param {object} root The container to work on.
 */
export function expectPlotButtons(root) {
  expect(within(root).getByRole('button', {name: 'Reset view'})).toBeInTheDocument()
  expect(within(root).getByRole('button', {name: 'Toggle fullscreen'})).toBeInTheDocument()
  expect(within(root).getByRole('button', {name: 'Capture image'})).toBeInTheDocument()
  expect(within(root).getByRole('button', {name: 'View data in the archive'})).toBeInTheDocument()
}
