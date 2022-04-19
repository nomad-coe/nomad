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
import { render, screen, startAPI, closeAPI } from '../conftest.spec'
import { navigateTo, browseRecursively } from '../archive/conftest.spec'
import EntryContext from './EntryContext'
import ArchiveEntryView from './ArchiveEntryView'

afterEach(() => closeAPI())

test('Browse archive reursively', async () => {
  await startAPI('tests.states.uploads.archive_browser_test', 'tests/data/uploads/archive_browser_test', 'test', 'password')
  const consoleLogSpy = jest.spyOn(console, 'log')
  const consoleErrorSpy = jest.spyOn(console, 'error')
  try {
    render(<EntryContext entryId={'1WGSYo1RrGFEIcM17Re4kjHC7k6p'}><ArchiveEntryView /></EntryContext>)
    expect(await screen.findByText('Entry')).toBeVisible()

    const path = ''
    const lane = await navigateTo(path)
    const laneIndex = path ? path.split('/').length : 0
    await browseRecursively(lane, laneIndex, join('*ArchiveBrowser*', path), consoleLogSpy, consoleErrorSpy)
  } finally {
    consoleLogSpy.mockRestore()
    consoleErrorSpy.mockRestore()
  }
}, 180000)
