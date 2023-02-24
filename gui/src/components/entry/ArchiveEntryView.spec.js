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
import { render, screen, within, startAPI, closeAPI, blockConsoleOutput, unblockConsoleOutput } from '../conftest.spec'
import { getLane, navigateTo, browseRecursively, pseudoRandomNumberGenerator } from '../archive/conftest.spec'
import { EntryContext } from './EntryContext'
import ArchiveEntryView from './ArchiveEntryView'
import { minutes } from '../../setupTests'
import { act } from 'react-dom/test-utils'
import {fireEvent} from '@testing-library/react'

let rand
/**
 * Lower values -> visits fewer locations.
 */
const visitProbabilityDecayFactor = 0.7

beforeEach(() => {
  blockConsoleOutput()
})
afterEach(() => {
  unblockConsoleOutput()
  closeAPI()
})

function archiveItemFilter(parentPath, itemKeys) {
  // The archive tree is very big and contains referential cycles, so we need to limit the crawling.
  // This method is used to make the selection.
  const segments = parentPath.split('/')
  if (segments.length === 1) {
    // Root - filter nothing
    return itemKeys
  }
  if (segments[segments.length - 2] === '_metainfo' || segments[segments.length - 2] === '_baseSectionDef@0') {
    // Never step deeper than one level into metainfo definitions, these are tested elsewhere
    // Never step deeper than one level into all inheriting secitons, the rest is tested seperately at the end
    return []
  }
  const rv = []
  const itemLists = {}
  for (const itemKey of itemKeys) {
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
  return rv.filter(() => rand() < visitProbabilityDecayFactor)
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

  render(<EntryContext entryId={entryId}><ArchiveEntryView /></EntryContext>)
  expect(await screen.findByText('Entry')).toBeVisible()

  if (withDefinition) {
    // Click the definitions checkbox
    await userEvent.click(screen.getByRoleAndText('checkbox', 'definitions'))
    expect(await within(getLane(0)).findByText('meta')).toBeVisible()
  }
  if (withAll) {
    // Click the metainfo definition
    await userEvent.click(screen.getByRoleAndText('checkbox', 'all defined'))
    expect(await within(getLane(0)).findByText('processing_logs')).toBeVisible()
  }
  await navigateTo(path)
  const laneIndex = path ? path.split('/').length : 0
  rand = pseudoRandomNumberGenerator()
  await browseRecursively(laneIndex, join(`*ArchiveBrowser ${name}*`, path), archiveItemFilter, filterKeyLength)
}, 20 * minutes) // NOTE!!! Do not increase this timeout! Rather, adjust the visitProbabilityDecayFactor

test('inheriting sections', async () => {
  await startAPI('tests.states.uploads.archive_browser_test', 'tests/data/uploads/archive_browser_test_inheriting_sectins', 'test', 'password')
  await act(async () => render(<EntryContext entryId={'Z0mBq-MtZ0B2IFveOhFCFJMPCZgO'}><ArchiveEntryView /></EntryContext>))
  expect(await screen.findByText('Entry')).toBeVisible()

  const path = 'data'
  const sectionName = '../uploads/archive_browser_test/raw/inheriting-schema.archive.yaml#/definitions/section_definitions/1'
  const sectionName2 = '../uploads/archive_browser_test/raw/inheriting-schema.archive.yaml#/definitions/section_definitions/2'
  await navigateTo(path)

  await userEvent.click(await screen.findByTestId('subsection:C1'))
  const selectLabel = await screen.findByText('Select a section')
  expect(selectLabel).toBeInTheDocument()

  const dropDown = await screen.findByTestId(`inheriting:SubSectionBase1`)
  expect(dropDown).toBeInTheDocument()
  const selectInput = within(dropDown).getByRole('textbox', { hidden: true })
  await waitFor(() => expect(selectInput.value).toEqual(sectionName))
  await fireEvent.change(selectInput, {target: {value: sectionName2}})
  await waitFor(() => expect(selectLabel).not.toBeInTheDocument())
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
  render(<EntryContext entryId={entryId}><ArchiveEntryView /></EntryContext>)
  expect(await screen.findByText('Entry')).toBeVisible()
  for (const path of paths) {
    await navigateTo(path)
  }
}, 5 * minutes)

test.each([
  [
    'author',
    'tests.states.entry.eln_properties',
    'tests/data/entry/eln_properties_writer',
    'test',
    'password'
  ], [
    'reviewer',
    'tests.states.entry.eln_properties',
    'tests/data/entry/eln_properties_reviewer',
    'ttester',
    'password'
  ]
])('ELN properties annotation as %s', async (name, state, snapshot, username, password) => {
  await startAPI(state, snapshot, username, password)
  await act(async () => render(<EntryContext entryId={'nwW7UiYFBvR9VQsKdehH9wd3Efws'}><ArchiveEntryView /></EntryContext>))
  expect(await screen.findByText('Entry')).toBeVisible()

  await navigateTo('data')

  await screen.getByText('Hall_measurement')
  const data = await screen.getByTestId('lane1:data')

  let compartments = await within(data).queryAllByTestId('compartment')
  const allVisibleQuantities = await within(data).queryAllByTestId("visible-quantity")
  const allVisibleSubsections = await within(data).queryAllByTestId("subsection")
  let allEditableQuantities = await within(data).queryAllByTestId("editable-quantity-editor")

  // test the visibility
  expect(compartments.length).toBe(4)
  expect(allVisibleQuantities.length).toBe(6)
  expect(allVisibleSubsections.length).toBe(3)

  // test the order
  await within(compartments[1]).findByText(/^quantities/i)
  await within(compartments[2]).findByText(/^sub section/i)
  await within(compartments[3]).findByText(/^referenced by/i)
  await within(allVisibleQuantities[0]).findByText(/^Crucible[ _]mass/i)
  await within(allVisibleQuantities[1]).findByText(/^Atmosphere/i)
  await within(allVisibleQuantities[2]).findByText(/^Crucible[ _]model/i)
  await within(allVisibleQuantities[3]).findByText(/^Brutto[ _]mass[ _]before/i)
  await within(allVisibleQuantities[4]).findByText(/^dose/i)
  await within(allVisibleQuantities[5]).findByText(/^net[ _]mass[ _]before/i)
  await within(allVisibleSubsections[0]).findByText(/^instrument2/i)
  await within(allVisibleSubsections[1]).findByText(/^instrument3/i)
  await within(allVisibleSubsections[2]).findByText(/^instrument1/i)

  // test the quantity permissions
  if (name === 'author') {
    expect(allEditableQuantities.length).toBe(3)
    await within(allEditableQuantities[0]).findByText(/^Crucible[ _]mass/i)
    await within(allEditableQuantities[1]).findByText(/^Atmosphere/i)
    await within(allEditableQuantities[2]).findByText(/^Brutto[ _]mass[ _]before/i)
  } else {
    expect(allEditableQuantities.length).toBe(0)
  }

  // test the subsection permissions
  await navigateTo('data/instrument2')
  const instrument2 = await screen.getByTestId('lane2:instrument2')
  compartments = await within(instrument2).queryAllByTestId('compartment')
  allEditableQuantities = await within(instrument2).queryAllByTestId("editable-quantity-editor")
  expect(compartments.length).toBe(4)
  expect(allEditableQuantities.length).toBe(0)

  // test the recursive subsection permissions
  await navigateTo('data/instrument2/internal1')
  const internal1 = await screen.getByTestId('lane3:internal1')
  compartments = await within(internal1).queryAllByTestId('compartment')
  allEditableQuantities = await within(internal1).queryAllByTestId("editable-quantity-editor")
  expect(compartments.length).toBe(3)
  expect(allEditableQuantities.length).toBe(0)
})
