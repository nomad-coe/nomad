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
import { fireEvent } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { waitFor } from '@testing-library/dom'
import { render, screen, within, startAPI, closeAPI } from '../conftest.spec'
import FileBrowser from './FileBrowser'
import { purgeTreePath, getLane, checkLanes, navigateTo, checkDirectoryLane, checkFileLane } from './conftest.spec'
import { minutes } from '../../setupTests'
import { createUploadUrl } from '../../utils'
import { apiBase } from '../../config'

const dirSpecialChars = 'dir special chars ~!?*\\()[]{}<>,.;:\'"`&@#$%=|'
const fileBrowserTree = {
  '': {cb: checkDirectoryLane},
  'deepdir': {cb: checkDirectoryLane},
  'deepdir/subdir1': {cb: checkDirectoryLane},
  'deepdir/subdir1/subdir2': {cb: checkDirectoryLane},
  'deepdir/subdir1/subdir2/dir special chars ~!?*\\()[]{}<>,.;:\'"`&@#$%=|': {cb: checkDirectoryLane},
  'deepdir/subdir1/subdir2/dir special chars ~!?*\\()[]{}<>,.;:\'"`&@#$%=|/f.txt': {cb: checkFileLane},
  'deepdir/subdir1/subdir2/dir special chars ~!?*\\()[]{}<>,.;:\'"`&@#$%=|/f.txt/preview': {cb: () => {}},
  'filetypes': {cb: checkDirectoryLane},
  'filetypes/image.png': {cb: checkFileLane},
  'filetypes/image.png/preview': {cb: () => {}},
  'filetypes/not an image.png': {cb: checkFileLane},
  'filetypes/not an image.png/preview': {cb: () => {}},
  'filetypes/doc.json': {cb: checkFileLane},
  'filetypes/doc.json/preview': {cb: () => {}},
  'test_entry': {cb: checkDirectoryLane},
  'test_entry/1.aux': {cb: checkFileLane},
  'test_entry/1.aux/preview': {cb: () => {}},
  'test_entry/2.aux': {cb: checkFileLane},
  'test_entry/2.aux/preview': {cb: () => {}},
  'test_entry/3.aux': {cb: checkFileLane},
  'test_entry/3.aux/preview': {cb: () => {}},
  'test_entry/4.aux': {cb: checkFileLane},
  'test_entry/4.aux/preview': {cb: () => {}},
  'test_entry/vasp.xml': {cb: checkFileLane, extra: {entryId: 'dft_bulk', parserName: 'parsers/vasp'}},
  'test_entry/vasp.xml/preview': {cb: () => {}}
}

afterEach(() => closeAPI())

async function testBrowseAround(editable) {
  // Create browser and check lanes
  const browserConfig = {
    rootTitle: 'Root Title',
    browserTree: fileBrowserTree,
    editable
  }
  render(<FileBrowser uploadUrl={createUploadUrl(apiBase, 'browser_test', '')} rootTitle={browserConfig.rootTitle}/>)
  await waitFor(() => {
    expect(screen.getByText('Root Title')).toBeVisible()
  })
  await checkLanes('', browserConfig)

  // Navigate around
  await navigateTo(`deepdir/subdir1/subdir2/${dirSpecialChars}/f.txt`, browserConfig)
  await within(getLane(5)).findByText(/preview/i) // auto-previewing .txt files with text viewer

  // Navigate to preview
  await navigateTo(`deepdir/subdir1/subdir2/${dirSpecialChars}/f.txt/preview`, browserConfig)
  await within(getLane(6)).findByText('content of f.txt') // auto-previewing .txt files with text viewer

  // Check file types

  // image: ok -> just check if there is an img tag
  await navigateTo('filetypes/image.png', browserConfig)
  await within(getLane(2)).findByText(/preview/i) // auto-previewing .txt files with text viewer

  await navigateTo('filetypes/image.png/preview', browserConfig)
  expect(within(getLane(3)).getByRole('img')).toBeVisible()

  // image: bad
  await navigateTo('filetypes/not an image.png', browserConfig)
  await within(getLane(2)).findByText(/preview/i) // auto-previewing .txt files with text viewer

  // Simulate a load error on the image (images are not actually loaded during testing, for efficiency reasons)
  await navigateTo('filetypes/not an image.png/preview', browserConfig)
  fireEvent.error(within(getLane(3)).getByRole('img'))
  await within(getLane(3)).findByText('Failed to open with image viewer. Bad file format?')
  expect(within(getLane(3)).getByButtonText('Open with text viewer')).toBeEnabled()
  await userEvent.click(within(getLane(3)).getByButtonText('Open with text viewer'))
  await within(getLane(3)).findByText('this is not an image!') // text content of the file

  // json
  await navigateTo('filetypes/doc.json', browserConfig)
  await within(getLane(2)).findByText(/preview/i)

  // json preview
  await navigateTo('filetypes/doc.json/preview', browserConfig)
  await within(getLane(3)).findByText('root')
  await within(getLane(3)).findByText(/"this is a json string value"/)

  // Check folder with an entry
  await navigateTo('test_entry/1.aux', browserConfig)
  await navigateTo('test_entry/vasp.xml', browserConfig)
  await within(getLane(2)).findByText(/preview/i)

  await navigateTo('test_entry/vasp.xml/preview', browserConfig)
  await within(getLane(3)).findByText(/GGA_X_PBE/) // This text should occur in the file preview
}

test.each([
  [
    'Browse unpublished as author',
    'tests.states.uploads.browser_test_unpublished',
    'tests/data/uploads/browser_test_unpublished_author',
    'test',
    'password',
    true
  ], [
    'Browse unpublished as reviewer',
    'tests.states.uploads.browser_test_unpublished',
    'tests/data/uploads/browser_test_unpublished_reviewer',
    'ttester',
    'password',
    false
  ], [
    'Browse published as author',
    'tests.states.uploads.browser_test_published',
    'tests/data/uploads/browser_test_published_author',
    'test',
    'password',
    false
  ]
])('Upload page: %s', async (name, state, snapshot, username, password, editable) => {
  await startAPI(state, snapshot, username, password)
  await testBrowseAround(editable)
}, 6 * minutes)

test('starting in entry dir', async () => {
  await startAPI('tests.states.uploads.browser_test_unpublished', 'tests/data/uploads/browser_test_entrydir', 'test', 'password')
  const entryDir = 'test_entry'
  // Extract subtree from the standard fileBrowserTree
  const fileBrowserTreeModified = {}
  for (let k in fileBrowserTree) {
    const v = fileBrowserTree[k]
    if (k === entryDir) {
      k = ''
    } else if (k.startsWith(entryDir + '/')) {
      k = k.substring(entryDir.length + 1)
    } else {
      continue
    }
    fileBrowserTreeModified[k] = v
  }

  const browserConfig = {
    rootTitle: 'Root Title',
    browserTree: fileBrowserTreeModified,
    editable: true
  }
  render(<FileBrowser uploadUrl={createUploadUrl(apiBase, 'browser_test', entryDir)} rootTitle={browserConfig.rootTitle}/>)
  await waitFor(() => {
    expect(screen.getByText('Root Title')).toBeVisible()
  })
  await checkLanes('', browserConfig)

  await navigateTo('1.aux', browserConfig)
  await navigateTo('vasp.xml', browserConfig)
  await within(getLane(1)).findByText(/preview/i) // This text should occur in the file preview

  await navigateTo('vasp.xml/preview', browserConfig)
  await within(getLane(2)).findByText(/GGA_X_PBE/) // This text should occur in the file preview
})

test('delete files', async () => {
  await startAPI('tests.states.uploads.browser_test_unpublished', 'tests/data/uploads/browser_test_delete_files', 'test', 'password')
  const fileBrowserTreeCopy = {...fileBrowserTree}
  const browserConfig = {
    rootTitle: 'Root Title',
    browserTree: fileBrowserTreeCopy,
    editable: true
  }
  render(<FileBrowser uploadUrl={createUploadUrl(apiBase, 'browser_test', '')} rootTitle={browserConfig.rootTitle}/>)
  await waitFor(() => {
    expect(screen.getByText('Root Title')).toBeVisible()
  })
  await checkLanes('', browserConfig)

  for (const fileName of ['vasp.xml', '1.aux']) {
    await navigateTo(`test_entry/${fileName}`, browserConfig)
    await userEvent.click(within(getLane(2)).getByButtonText('delete this file'))
    await screen.findByText(/Really delete the file/)
    expect(screen.queryAllByText(fileName).length).toBeGreaterThan(1)
    await userEvent.click(screen.getByButtonText('OK'))
    await waitFor(() => {
      expect(screen.queryAllByText(fileName).length).toEqual(0)
    })
    purgeTreePath(fileBrowserTreeCopy, `test_entry/${fileName}`)
    checkLanes('test_entry', browserConfig)
  }
})

test('delete folder', async () => {
  await startAPI('tests.states.uploads.browser_test_unpublished', 'tests/data/uploads/browser_test_delete_folders', 'test', 'password')
  const fileBrowserTreeCopy = {...fileBrowserTree}
  const browserConfig = {
    rootTitle: 'Root Title',
    browserTree: fileBrowserTreeCopy,
    editable: true
  }
  render(<FileBrowser uploadUrl={createUploadUrl(apiBase, 'browser_test', '')} rootTitle={browserConfig.rootTitle}/>)
  await waitFor(() => {
    expect(screen.getByText('Root Title')).toBeVisible()
  })
  await checkLanes('', browserConfig)

  for (const path of ['test_entry', 'deepdir/subdir1/subdir2']) {
    const segments = path.split('/')
    const folderName = segments[segments.length - 1]
    const parentPath = (segments.slice(0, segments.length - 1)).join('/')
    // Navigate to path
    const lane = await navigateTo(path, browserConfig)
    await userEvent.click(within(lane).getByButtonText('delete this folder'))
    await screen.findByText(/Really delete the directory/)
    expect(screen.queryAllByText(folderName).length).toBeGreaterThan(1)
    await userEvent.click(screen.getByButtonText('OK'))
    await waitFor(() => {
      expect(screen.queryAllByText(folderName).length).toEqual(0)
    })
    purgeTreePath(fileBrowserTreeCopy, path)
    checkLanes(parentPath, browserConfig)
  }
})
