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
import userEvent from '@testing-library/user-event'
import { render, screen, within, startAPI, closeAPI, blockConsoleOutput, unblockConsoleOutput, waitForGUI } from '../conftest.spec'
import { getLane, navigateTo, browseRecursively } from '../archive/conftest.spec'
import EntryPageContext from './EntryPageContext'
import ArchiveEntryView from './ArchiveEntryView'
import { minutes } from '../../setupTests'

beforeEach(() => {
  blockConsoleOutput()
})
afterEach(() => {
  unblockConsoleOutput()
  closeAPI()
})

function archiveItemFilter(parentPath, items) {
  // The archive tree is very big and contains referential cycles, so we need to limit the crawling.
  // This method is used to make the selection.
  const segments = parentPath.split('/')
  if (segments.length === 1) {
    // Root - filter nothing
    return Object.keys(items)
  }
  if (segments[segments.length - 2] === '_metainfo' || segments[segments.length - 2] === '_baseSectionDef@0') {
    // Never step deeper than one level into metainfo definitions, these are tested elsewhere
    // Never step deeper than one level into all inheriting secitons, the rest is tested seperately at the end
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
  ['vasp', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', '', false, false, 1],
  ['vasp with definitions', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', '', true, false, 2],
  ['vasp with all', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', '', false, true, 1],
  ['Sample', '6x6VbK15sTesOX3wHPHp2zTiEo8Y', '', false, false, 2],
  ['Sample with definitions', '6x6VbK15sTesOX3wHPHp2zTiEo8Y', '', true, false, 2],
  ['Chemical', '8mYi1m21s3k2LdjIUt0EnLyQ-_-c', '', false, false, 2],
  ['Schema', 'EdNTmfG89xq927d5Zcq-hMlLWaa2', '', false, false, 2],
  ['MySection-inter-entry', 'QoV6GTCH0lT4ArliEfO6ytXa9ihg', '', false, false, 2],
  ['MySection-intra-entry', 'oqOTEw8ZzGt4Z763KE2IUv9iY1M5', '', false, false, 2],
  ['With-inheriting-sections', 'Z0mBq-MtZ0B2IFveOhFCFJMPCZgO', '', false, false, 2],
  ['With-inheriting-sections with definitions', 'Z0mBq-MtZ0B2IFveOhFCFJMPCZgO', '', true, false, 2]
])('Browse archive recursively: %s', async (name, entryId, path, withDefinition, withAll, filterKeyLength) => {
  await startAPI(
    'tests.states.uploads.archive_browser_test',
    'tests/data/uploads/archive_browser_test_' + name.replace(/ /g, '_'),
    'test', 'password')

  render(<EntryPageContext entryId={entryId}><ArchiveEntryView /></EntryPageContext>)
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
  await browseRecursively(lane, laneIndex, join(`*ArchiveBrowser ${name}*`, path), archiveItemFilter, filterKeyLength)
}, 20 * minutes)

test('inheriting sections', async () => {
  await startAPI('tests.states.uploads.archive_browser_test', 'tests/data/uploads/archive_browser_test_inheriting_sectins', 'test', 'password')

  render(<EntryPageContext entryId={'Z0mBq-MtZ0B2IFveOhFCFJMPCZgO'}><ArchiveEntryView /></EntryPageContext>)
  expect(await screen.findByText('Entry')).toBeVisible()

  const path = 'data'
  const sectionName = '../uploads/archive_browser_test/raw/inheriting-schema.archive.yaml#definitions/section_definitions/1'
  await navigateTo(path)

  userEvent.click(await screen.findByTestId('subsection:C1'))
  expect(await screen.findByText('Select an m_def from the list')).toBeInTheDocument()

  const dropDown = await screen.findByTestId(`inheriting:SubSectionBase1`)
  expect(dropDown).toBeInTheDocument()
  const selectInput = within(dropDown).getByRole('textbox', { hidden: true })

  await waitForGUI()
  await waitFor(() => expect(selectInput.value).toEqual(`${sectionName}`))
})

test.each([
  ['vasp', '1WGSYo1RrGFEIcM17Re4kjHC7k6p', [
    'metadata/last_processing_time',
    'results/material/elements',
    'run/calculation/stress/total']]
])('specific paths %s', async (name, entryId, paths) => {
  // Test some paths that we expect to exist
  await startAPI(
    'tests.states.uploads.archive_browser_test',
    'tests/data/uploads/archive_browser_test_paths_' + name.replace(/ /g, '_'),
    'test', 'password')
  render(<EntryPageContext entryId={entryId}><ArchiveEntryView /></EntryPageContext>)
  expect(await screen.findByText('Entry')).toBeVisible()
  for (const path of paths) {
    await navigateTo(path)
  }
}, 5 * minutes)
