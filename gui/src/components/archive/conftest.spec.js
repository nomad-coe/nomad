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
import { waitFor } from '@testing-library/dom'
import { screen, within } from '../conftest.spec'
import userEvent from '@testing-library/user-event'

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
      lanePath += (lanePath ? '/' : '') + segment
      lanePaths.push(lanePath)
    })
  }
  for (let laneIndex = 0; laneIndex < lanePaths.length; laneIndex++) {
    const lanePath = lanePaths[laneIndex]
    const { cb, extra } = browserTree[lanePath]
    expect(cb).toBeDefined()
    const lastSegment = lanePath.split('/').pop()
    const lane = getLane(laneIndex, lastSegment)
    const args = {lane, laneIndex, lanePath, lastSegment, ...browserConfig, ...extra}
    await cb(args)
  }
  expect(getLane(lanePaths.length)).toBeNull() // otherwise we have too many lanes
}

/**
 * Navigates to the specified path by firing a click event on the last item in the path.
 * The parent lane (i.e. the lane containing the item to click) must already be open in the
 * browser. After firing the click event, waits for the rendering to complete, and check
 * the lanes.
 */
export async function navigateAndCheck(path, browserConfig) {
  const segments = path.split('/')
  const item = segments[segments.length - 1]
  const parentLaneIndex = segments.length - 1
  userEvent.click(within(getLane(parentLaneIndex)).getByText(item))
  await waitFor(() => {
    expect(getLane(parentLaneIndex + 1, item)).not.toBeNull()
    expect(getLane(parentLaneIndex + 2)).toBeNull()
  })
  await checkLanes(path, browserConfig)
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
export async function checkFileLane(
  {lane, lastSegment, entryId, parserName, editable}) {
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
