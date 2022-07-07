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
import { join } from 'path'
import { waitFor } from '@testing-library/dom'
import { render, screen, blockConsoleOutput, unblockConsoleOutput } from '../conftest.spec'
import { navigateTo, browseRecursively } from './conftest.spec'
import MetainfoBrowser from './MetainfoBrowser'
import { minutes } from '../../setupTests'

function metainfoItemFilter(parentPath, items) {
  // The metainfo tree is very big, so we need to limit the crawling. This method is used
  // to make the selection.
  if (parentPath.split('/').length === 1) {
    // Root - filter nothing
    return Object.keys(items)
  }
  const segments = parentPath.split('/')
  if (segments[segments.length - 1].includes('_inheritingSectionDef@') || segments[segments.length - 1].includes('_baseSectionDef@')) {
    // Never step into allInheritingSections
    return []
  }
  const rv = []
  const encounteredPrefixes = {}
  let countKeysWithoutPrefixes = 0
  for (const itemKey of Object.keys(items)) {
    let include = false
    const parts = itemKey.split(/[.@]/)
    const prefixes = parts.slice(0, parts.length - 1)
    if (prefixes.length === 0) {
      countKeysWithoutPrefixes++
      include = countKeysWithoutPrefixes < 3
    } else {
      // We have some number of prefixes. Include only if at least one of the prefixes is new.
      for (const prefix of prefixes) {
        if (!encounteredPrefixes[prefix]) {
          include = true
          encounteredPrefixes[prefix] = true
        }
      }
    }
    if (include) {
      rv.push(itemKey)
    }
  }
  return rv
}

beforeEach(() => blockConsoleOutput())
afterEach(() => unblockConsoleOutput())

test('Browse metainfo reursively', async () => {
  render(<MetainfoBrowser />)
  await waitFor(() => {
    expect(screen.getByText(/archive root section/i)).toBeVisible()
  })

  const path = ''
  const lane = await navigateTo(path)
  const laneIndex = path ? path.split('/').length : 0
  await browseRecursively(lane, laneIndex, join('*MetaInfoBrowser*', path), metainfoItemFilter, 2)
}, 30 * minutes)
