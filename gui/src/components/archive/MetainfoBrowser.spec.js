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
import { navigateTo, browseRecursively, pseudoRandomNumberGenerator } from './conftest.spec'
import MetainfoBrowser from './MetainfoBrowser'
import { minutes } from '../../setupTests'

const rand = pseudoRandomNumberGenerator()

/**
 * Lower values -> visits fewer locations.
 */
const visitProbabilityDecayFactor = 0.1

function metainfoItemFilter(parentPath, itemKeys) {
  // The metainfo tree is very big, so we need to limit the crawling. This method is used
  // to make the selection.
  const parentSegments = parentPath.split('/').length - 1
  if (parentSegments === 0) {
    // Root - filter nothing
    return itemKeys
  }
  const rv = []
  const categoryCache = {}
  for (const itemKey of itemKeys) {
    // Compute a "category" for the item, which is a string based on certain properties of
    // the itemKey. When selecting items to visit, we try to cover as many categories as possible
    let category = itemKey.startsWith('_') ? '1' : '0'
    category += itemKey.includes('_') ? '1' : '0'
    category += itemKey.includes('@') ? '1' : '0'
    category += itemKey.match(/A-Z/) ? '1' : '0'
    category += itemKey.startsWith('section_definitions') ? '1' : '0'
    category += itemKey.startsWith('category_definitions') ? '1' : '0'
    let cache = categoryCache[category]
    if (!cache) {
      cache = []
      categoryCache[category] = cache
    }
    cache.push(itemKey)
  }
  const visitProbability = visitProbabilityDecayFactor ** parentSegments
  for (const category of Object.keys(categoryCache).sort()) {
    if (rand() < visitProbability) {
      // Include one of the itemKeys in this category, selected at random
      const categoryItems = categoryCache[category]
      rv.push(categoryItems[Math.floor(rand() * categoryItems.length)])
    }
  }
  return rv.sort()
}

beforeEach(() => blockConsoleOutput())
afterEach(() => unblockConsoleOutput())

test('Browse metainfo pseudorandomly', async () => {
  render(<MetainfoBrowser />)
  await waitFor(() => {
    expect(screen.getByText(/archive root section/i)).toBeVisible()
  })

  const path = ''
  await navigateTo(path)
  const laneIndex = path ? path.split('/').length : 0
  const {count} = await browseRecursively(laneIndex, join('*MetaInfoBrowser*', path), metainfoItemFilter, 2)

  // Currently we do not test that the visited items are the same using the hash
  // returned by browseRecursively. This is because any change in the metainfo
  // schema may modify the traversal order, thus changing the hash. If a fixed
  // traversal would be required, the list of traversed items should be given
  // explicitly, but that then can break if any of the items is removed/renamed.

  // Check that the tested number of paths is enough, but also not too high. Adjust
  // visitProbabilityDecayFactor if the number is not in this range.
  expect(count).toBeGreaterThan(100)
  expect(count).toBeLessThan(700)
}, 20 * minutes) // NOTE!!! Do not increase this timeout! Rather, adjust the visitProbabilityDecayFactor
