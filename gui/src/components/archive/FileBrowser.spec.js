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
import { render, screen, wait } from '../../testSetup'
import FileBrowser from './FileBrowser'
import { useApi } from '../api'
import { waitFor, within } from '@testing-library/dom'
import userEvent from '@testing-library/user-event'

jest.mock('../api')

const dirSpecialChars = 'dir special chars ~!?*\\()[]{}<>,.;:\'"`&@#$%=|'
const directoryTree = {
  '': [
    {name: 'rootdir1', is_file: false, size: 456},
    {name: 'rootdir2', is_file: false, size: 456},
    {name: 'rootfile1', is_file: true, size: 123}],
  'rootdir1': [
    {name: 'subdir1', is_file: false, size: 456},
    {name: 'file1.json', is_file: true, size: 123, entry_id: 'entry_1', parser_name: 'parsers/vasp'}],
  'rootdir2': [
    {name: dirSpecialChars, is_file: false, size: 123}],
  'rootdir2/dir special chars ~!?*\\()[]{}<>,.;:\'"`&@#$%=|': [
    {name: 'file2', is_file: true, size: 123}]
}

beforeAll(() => {
  // API mock init
  useApi.mockReturnValue({
    api: {
      get: (url) => {
        const baseUrl = url.split('?')[0]
        const decodedBaseUrl = baseUrl.split('/').map(segment => decodeURIComponent(segment)).join('/')
        const path = decodedBaseUrl.split('/').slice(4).join('/')
        const segments = path.split('/')
        const dirContent = directoryTree[path]
        if (dirContent === undefined) {
          return wait({detail: 'Not found. Invalid path?'})
        }
        const response = {
          path: path,
          access: 'unpublished',
          directory_metadata: {
            name: segments[segments.length - 1],
            size: 456,
            content: dirContent.map(e => { return {...e} })
          },
          pagination: {
            page_size: 500,
            page: 1,
            total: dirContent.length
          }
        }
        return wait(response)
      }
    }
  })
})

afterAll(() => {
  // API mock cleanup
  jest.unmock('../api')
})

function checkLanes(path, rootPath, rootTitle) {
  // Checks the content of the lanes
  let lanePath = rootPath
  const lanePaths = [lanePath]
  path.split('/').forEach(segment => {
    if (segment) {
      lanePath += (lanePath ? '/' : '') + segment
      lanePaths.push(lanePath)
    }
  })
  lanePaths.forEach((lanePath, i) => {
    const segments = lanePath.split('/')
    const lastSegment = segments[segments.length - 1]
    const lane = screen.getByTestId(`lane${i}`)
    const expectedContent = directoryTree[lanePath]
    if (expectedContent === undefined) {
      // File lane
      expect(i).toBeGreaterThan(0)
      const directoryData = directoryTree[lanePaths[i - 1]]
      expect(directoryData).toBeDefined()
      const fileData = directoryData.filter(e => e.name === lastSegment)[0]
      expect(fileData).toBeDefined()
      expect(within(lane).getByText(lastSegment)).toBeVisible() // Lane title = file name
      if (fileData.entry_id) {
        expect(within(lane).getByText(fileData.entry_id)).toBeVisible()
        expect(within(lane).getByText(fileData.parser_name.replace('parser/', ''))).toBeVisible()
      }
    } else {
      // Directory lane
      expect(within(lane).getByText(i === 0 ? rootTitle : lastSegment)).toBeVisible()
      expectedContent.forEach(e => {
        expect(within(lane).getByText(e.name)).toBeVisible()
      })
    }
  })
  expect(screen.queryByTestId(`lane${lanePaths.length}`)).not.toBeInTheDocument()
}

test('render browser and browse around', async () => {
  render(<FileBrowser uploadId="upload_id_1" path="" rootTitle="Root Title"/>)
  await waitFor(() => {
    expect(screen.getByText('Root Title')).toBeVisible()
  })
  checkLanes('', '', 'Root Title')
  // Select item: rootdir1
  userEvent.click(screen.getByText('rootdir1'))
  await waitFor(() => {
    expect(screen.queryByTestId('lane1')).toBeInTheDocument()
  })
  checkLanes('rootdir1', '', 'Root Title')
  // Select item: rootdir1/file1.json
  userEvent.click(screen.getByText('file1.json'))
  await waitFor(() => {
    expect(screen.queryByTestId('lane2')).toBeInTheDocument()
  })
  checkLanes('rootdir1/file1.json', '', 'Root Title')
  // Select item: rootdir2
  userEvent.click(screen.getByText('rootdir2'))
  await waitFor(() => {
    expect(within(screen.queryByTestId('lane1')).getByText('rootdir2')).toBeInTheDocument()
  })
  checkLanes('rootdir2', '', 'Root Title')
  // Select item: rootdir2/dirSpecialChars
  userEvent.click(screen.getByText(dirSpecialChars))
  await waitFor(() => {
    expect(screen.queryByTestId('lane2')).toBeInTheDocument()
  })
  checkLanes(`rootdir2/${dirSpecialChars}`, '', 'Root Title')
  // Select item: rootdir2/dirSpecialChars/file2
  userEvent.click(screen.getByText('file2'))
  await waitFor(() => {
    expect(screen.queryByTestId('lane3')).toBeInTheDocument()
  })
  checkLanes(`rootdir2/${dirSpecialChars}/file2`, '', 'Root Title')
})
