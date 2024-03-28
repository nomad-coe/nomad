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
import { join, basename } from 'path'
import { waitFor } from '@testing-library/dom'
import { screen, within, expectNoConsoleOutput } from '../conftest.spec'
import userEvent from '@testing-library/user-event'
import { Item, Content, Compartment, Title, laneErrorBoundryMessage, Adaptor } from './Browser'
import { isWaitingForUpdateTestId } from '../../utils'
import React from 'react'
import PropTypes from 'prop-types'

const crypto = require('crypto')

/*****************************************************************************************
 * Utilities for testing browser functionality.
 *
 * Key concepts:
 *  browserConfig
 *    An object storing information about how a Browser component was instantiated and configured,
 *    plus a *browserTree* attriute defining what should be shown in the browser for different
 *    paths (depending on where the user navigates).
 *  browserTree
 *    An object with keys and values, where the keys are browseable paths ('' = the browser root lane)
 *    and the values are objects of the form:
 *      {cb: <lane check function>, extra: <object with additional key-value pairs>}
 *  lane check function
 *    A function which performs the standard checks that a particular lane is rendered correctly.
 *    These functions are referred in the *browserTree* (the *cb* attribute), and the checkLanes
 *    utility function uses them to check the lanes that should be visible for a certain path.
 *    The function is invoked with the object
 *      {lane, laneIndex, lanePath, lastSegment, ...browserConfig, ...extra}
 *    as input. Some lane check functions are defined in this file, as they may be used in
 *    tests defined in different files.
 *****************************************************************************************/

/**
 * Utility for calculating which items should exist in the given path in the specified
 * browserTree. Returns these as an array of strings.
 */
export function itemsInTreePath(browserTree, path) {
  const rv = []
  const pathPrefix = path ? path + '/' : ''
  for (const treePath in browserTree) {
    if (treePath && treePath.startsWith(pathPrefix)) {
      const item = treePath.substring(pathPrefix.length).split('/')[0]
      if (!rv.includes(item)) {
        rv.push(item)
      }
    }
  }
  return rv.sort()
}

/**
 * Deletes everything at or under a specific path in the provided browserTree
 */
export function purgeTreePath(browserTree, path) {
  const keysToDelete = []
  for (const key in browserTree) {
    if (!path || key === path || key.startsWith(path + '/')) {
      keysToDelete.push(key)
    }
  }
  keysToDelete.forEach(key => delete browserTree[key])
}

/**
 * Gets a lane object from the screen by its index. If a laneKey is specified, we also check
 * that the lane has this key. If no matching lane is found, we return null.
 */
export function getLane(laneIndex, laneKey) {
  try {
    if (laneKey) {
      laneKey = laneKey.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') // Escape for use in regexp
    }
    const re = laneKey ? new RegExp(`^lane${laneIndex}:${laneKey}`) : new RegExp(`^lane${laneIndex}:.*`)
    return screen.getByTestId(re)
  } catch {
    return null
  }
}

/**
 * Extracts the lane key from a lane object (using its data-testid attribute)
 */
export function getLaneKey(lane) {
  const testId = lane.attributes['data-testid'].value
  return testId.substring(testId.search(':') + 1)
}

/**
 * Extracts the item key from an item object (using its data-testid attribute)
 */
export function getItemKey(item) {
  const testId = item.attributes['data-testid'].value
  return testId.substring(5)
}

/**
 * Checks that the browser is correctly rendered, provided the browser's path.
 * Each lane corresponding to a segment in the path is checked by calling the cb function
 * defined in the browserTree (which should be an attribute of browserConfig)
 */
export async function checkLanes(path, browserConfig) {
  const { browserTree } = browserConfig
  const lanePaths = ['']
  if (path) {
    let lanePath = ''
    path.split('/').forEach(segment => {
      lanePath = join(lanePath, segment)
      lanePaths.push(lanePath)
    })
  }
  for (let laneIndex = 0; laneIndex < lanePaths.length; laneIndex++) {
    const lanePath = lanePaths[laneIndex]
    const { cb, extra } = browserTree[lanePath]
    expect(cb).toBeDefined()
    const lastSegment = basename(lanePath)
    const lane = getLane(laneIndex, lastSegment)
    const args = {lane, laneIndex, lanePath, lastSegment, ...browserConfig, ...extra}
    await cb(args)
  }
  expect(getLane(lanePaths.length)).toBeNull() // otherwise we have too many lanes
}

/**
 * Selects (clicks) on an item in the specified lane and waits for the new lane to render.
 * Note, the item must not already be selected. Throws an exception if the rendering fails.
 * Returns the new lane when done. If you already have the item object, you can pass it
 * as the last argument (for optimization purposes), otherwise it will be fetched from the
 * itemKey.
 */
export async function selectItemAndWaitForRender(laneIndex, itemKey) {
  const lane = getLane(laneIndex)
  const item = within(lane).getByTestId(`item:${itemKey}`)
  await userEvent.click(item)
  await waitFor(() => {
    const nextLane = getLane(laneIndex + 1, itemKey)
    expect(nextLane).not.toBeNull()
    expect(getLane(laneIndex + 2)).toBeNull()
    expect(within(nextLane).queryAllByTestId(isWaitingForUpdateTestId).length).toBe(0)
  })
  const nextLane = getLane(laneIndex + 1, itemKey)
  expect(within(nextLane).queryByText(laneErrorBoundryMessage)).toBeNull()
  return nextLane
}

/**
 * Navigates to the specified path by clicking the appropriate items in the lanes. The method
 * figures out itself, by inspecting the existing lanes, which clicks are needed. If a browserConfig
 * with a browserTree attribute is specified, we also call checkLanes after each click.
 * Returns the final lane.
 */
export async function navigateTo(path, browserConfig) {
  if (!path) {
    return getLane(0)
  }
  const segments = path.split('/')
  let subpath = ''
  let laneIndex = 0
  let lane
  let nextLane = getLane(1)
  for (const segment of segments) {
    subpath = join(subpath, segment)
    const selectedItemKey = nextLane ? getLaneKey(nextLane) : null
    if (selectedItemKey !== segment) {
      // Need to select a (different) item in this lane
      lane = await selectItemAndWaitForRender(laneIndex, segment)
      if (browserConfig?.browserTree) {
        await checkLanes(subpath, browserConfig)
      }
    } else {
      lane = nextLane
    }
    laneIndex += 1
    nextLane = getLane(laneIndex + 1)
  }
  return lane
}

/**
 * Browses recursively, navigating to all items that pass the provided itemFilter. If no
 * itemFilter is provided, all items will be visited. This method provides an easy way to
 * verify that the browser renders all (or at least a lot of) paths correctly.
 */
export async function browseRecursively(laneIndex, path, itemFilter, filterKeyLength = 2, filterMemory = null) {
  const lane = getLane(laneIndex)
  let count = 0
  const hash = crypto.createHash('sha512')

  if (filterMemory === null) {
    filterMemory = new Set()
  }
  // Click on all discovered item-lists to open them
  for (const itemList of within(lane).queryAllByRole('item-list')) {
    const label = itemList.textContent
    try {
      await userEvent.click(itemList)
      await within(lane).findByTestId(`item-list:${label}`)
      expectNoConsoleOutput()
    } catch (error) {
      process.stdout.write(`ERROR expanding item list: ${path}/${label}\n`)
      throw error
    }
  }
  // Click on all collapsed compartments to uncollapse them
  const collapsedChips = within(lane).queryAllByTestId(/^collapsed:/)
  for (let idx = 0; idx < collapsedChips.length; idx++) {
    try {
      await userEvent.click(collapsedChips[idx])
      await waitFor(() => {
        expect(within(lane).queryAllByTestId(/^collapsed:/).length).toBe(collapsedChips.length - idx - 1)
        expect(within(lane).queryAllByTestId(isWaitingForUpdateTestId).length).toBe(0)
      })
    } catch (error) {
      process.stdout.write(`ERROR uncollapsing compartment ${idx}`)
      throw error
    }
  }
  let itemKeys = []
  for (const item of within(lane).queryAllByTestId(/^item:/)) {
    const itemKey = getItemKey(item)
    itemKeys.push(itemKey)
    hash.update(itemKey)
  }
  if (itemFilter) {
    itemKeys = itemFilter(path, itemKeys)
  }
  for (const itemKey of itemKeys) {
    const itemPath = join(path, itemKey)
    const segments = itemPath.split('/')
    const filterKey = segments.slice(segments.length - filterKeyLength).join('/')
    if (!filterMemory.has(filterKey)) {
      filterMemory.add(filterKey)
      const nextPath = `${path}/${itemKey}`
      // Uncomment line below if you want to show the paths visited.
      // process.stdout.write(`next path: ${nextPath}\n`)
      // process.stdout.write(`mem: ${process.memoryUsage().heapTotal / 1e9}\n`)
      try {
        await selectItemAndWaitForRender(laneIndex, itemKey)
        expectNoConsoleOutput()
      } catch (error) {
        process.stdout.write(`ERROR encountered when browsing to: ${nextPath}\n`)
        throw error
      }
      // new lane rendered successfully
      count++
      const rv = await browseRecursively(laneIndex + 1, nextPath, itemFilter, filterKeyLength, filterMemory)
      count += rv.count
      hash.update(rv.hash)
    }
  }
  return {count, hash: hash.digest('base64').slice(0, 28)}
}

/*****************************************************************************************
 * Lane check functions
 *****************************************************************************************/

/**
 * Lane check function for directory lanes
 */
export async function checkDirectoryLane({lane, laneIndex, lanePath, lastSegment, browserTree, rootTitle, editable}) {
  expect(within(lane).getByText(laneIndex === 0 ? rootTitle : lastSegment)).toBeVisible() // Lane title

  itemsInTreePath(browserTree, lanePath).forEach(item => {
    expect(within(lane).getByText(item)).toBeVisible()
  })
  // Buttons
  expect(within(lane).getByButtonText('download this folder')).toBeEnabled()
  expect(within(lane).getByButtonText('reload directory contents')).toBeEnabled()
  for (const buttonTitle of ['upload to this folder (click or drop files)', 'create new folder', 'delete this folder']) {
    if (editable) {
      expect(within(lane).getByButtonText(buttonTitle)).toBeEnabled()
    } else {
      expect(within(lane).queryByButtonText(buttonTitle)).toBeNull()
    }
  }
}

/**
 * Lane check function for file preview lanes
 */
export async function checkFileLane({lane, lastSegment, entryId, parserName, editable}) {
  expect(within(lane).getByText(lastSegment)).toBeVisible() // Lane title
  if (entryId) {
    expect(within(lane).getByText(entryId)).toBeVisible()
  }
  if (parserName) {
    expect(within(lane).getByText(parserName)).toBeVisible()
  }
  // Buttons
  expect(within(lane).getByButtonText('download this file')).toBeEnabled()
  if (editable) {
    expect(within(lane).getByButtonText('delete this file')).toBeEnabled()
  } else {
    expect(within(lane).queryByButtonText('delete this file')).toBeNull()
  }
}

/*****************************************************************************************
 * Misc
 *****************************************************************************************/

/**
 * A simple pseudo-random number generator
 */
export function pseudoRandomNumberGenerator(seed = 7) {
  return () => {
    seed += 0x6D2B79F5
    let t = seed
    t = Math.imul(t ^ t >>> 15, t | 1)
    t ^= t + Math.imul(t ^ t >>> 7, t | 61)
    return ((t ^ t >>> 14) >>> 0) / 4294967296
  }
}

function checkTestLane({lane, laneIndex, lanePath, lastSegment, browserTree, rootTitle}) {
  expect(within(lane).getByText(laneIndex === 0 ? rootTitle : lastSegment)).toBeVisible() // Lane title

  itemsInTreePath(browserTree, lanePath).forEach(item => {
    expect(within(lane).getByText(item)).toBeVisible()
  })
}

export const browserTree = {
  '': {cb: checkTestLane},
  'dir1': {cb: checkTestLane},
  'dir1/success': {cb: checkTestLane},
  'dir1/fail': {cb: checkTestLane}
}

export class TestAdaptor extends Adaptor {
  constructor(path, title) {
    super()
    this.path = path
    this.title = title
    this.parsedObjUrl = {entryId: 'entryID1'}
  }

  itemAdaptor(key) {
    return new TestAdaptor(join(this.path, key))
  }

  render() {
    return <TestContent path={this.path} title={this.title} />
  }
}

export function TestContent({path, title}) {
  if (path.endsWith('/fail')) {
    throw new Error('mocked render error')
  }
  return <Content key={path}>
    <Title title={title || basename(path)} label="test" />
    <Compartment>
      {
        itemsInTreePath(browserTree, path).map(itemKey => (
          <Item itemKey={itemKey} key={join(path, itemKey)}>
            {itemKey}
          </Item>
        ))
      }
    </Compartment>
  </Content>
}
TestContent.propTypes = {
  path: PropTypes.string.isRequired,
  title: PropTypes.string
}
