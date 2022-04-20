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
import userEvent from '@testing-library/user-event'
import { render, screen, within, startAPI, closeAPI } from '../conftest.spec'
import { getLane, navigateTo, browseRecursively } from '../archive/conftest.spec'
import EntryContext from './EntryContext'
import ArchiveEntryView from './ArchiveEntryView'
import { minutes } from '../../setupTests'

afterEach(() => closeAPI())

function archiveItemFilter(parentPath, items) {
  // The archive tree is very big and contains referential cycles, so we need to limit the crawling.
  // This method is used to make the selection.
  const segments = parentPath.split('/')
  if (segments.length === 1) {
    // Root - filter nothing
    return Object.keys(items)
  }
  if (segments[segments.length - 2] === '_metainfo') {
    // Never step deeper than one level into metainfo definitions, these are tested elsewhere
    return []
  }
  const rv = []
  const itemLists = {}
  for (const itemKey of Object.keys(items)) {
    const parts = itemKey.split(':')
    if (parts.length === 2) {
      const [label, index] = parts
      if (itemLists[label]) {
        itemLists[label].push(index)
      } else {
        itemLists[label] = [index]
      }
    } else {
      rv.push(itemKey)
    }
  }
  // For item lists, we want to browse only the first and the last index
  for (const [label, indices] of Object.entries(itemLists)) {
    rv.push(`${label}:${indices[0]}`)
    rv.push(`${label}:${indices[indices.length - 1]}`)
  }
  return rv
}

test.each([
  ['normal', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', '', false, false, 1],
  ['with definitions', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', '', true, false, 2],
  ['with all', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', '', false, true, 1]
])('Browse archive recursively: %s', async (name, entryId, path, withDefinition, withAll, filterKeyLength) => {
  await startAPI('tests.states.uploads.archive_browser_test', 'tests/data/uploads/archive_browser_test', 'test', 'password')
  const consoleLogSpy = jest.spyOn(console, 'log')
  const consoleErrorSpy = jest.spyOn(console, 'error')
  try {
    render(<EntryContext entryId={entryId}><ArchiveEntryView /></EntryContext>)
    expect(await screen.findByText('Entry')).toBeVisible()

    if (withDefinition) {
      // Click the definitions checkbox
      userEvent.click(screen.getByRoleAndText('checkbox', 'definitions'))
      expect(await within(getLane(0)).findByText('meta')).toBeVisible()
    }
    if (withAll) {
      // Click the metainfo definition
      userEvent.click(screen.getByRoleAndText('checkbox', 'all defined'))
      expect(await within(getLane(0)).findByText('processing_logs')).toBeVisible()
    }
    const lane = await navigateTo(path)
    const laneIndex = path ? path.split('/').length : 0
    await browseRecursively(lane, laneIndex, join(`*ArchiveBrowser ${name}*`, path), consoleLogSpy, consoleErrorSpy, archiveItemFilter, filterKeyLength)
  } finally {
    consoleLogSpy.mockRestore()
    consoleErrorSpy.mockRestore()
  }
}, 12 * minutes)
