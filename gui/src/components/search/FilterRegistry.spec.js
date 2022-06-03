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

import { FilterRegistry } from './FilterRegistry'

test('test abbreviations', async () => {
  const registry = new FilterRegistry()
  registry.register('a.b', {})

  // No overlap in abbreviations
  expect(registry.abbreviations['a.b']).toBe('b')
  expect(registry.fullnames['b']).toBe('a.b')

  // After overlap in abbreviations
  registry.register('a.b.b', {})
  expect(registry.abbreviations['a.b']).toBe('a.b')
  expect(registry.abbreviations['a.b.b']).toBe('a.b.b')
  expect(registry.fullnames['b']).toBe(undefined)
})
