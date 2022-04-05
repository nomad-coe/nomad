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

/* eslint-disable no-undef */

import React from 'react'
import { get, isPlainObject } from 'lodash'
import { rest } from 'msw'
import PropTypes from 'prop-types'
import { RecoilRoot } from 'recoil'
import {
  render,
  screen,
  within,
  buildQueries,
  queryAllByText
} from '@testing-library/react'
import { server } from '../setupTests'
import { getArchive } from '../../tests/DFTBulk'
import { Router, MemoryRouter } from 'react-router-dom'
import { KeycloakProvider } from 'react-keycloak'
import { createBrowserHistory } from 'history'
import { APIProvider } from './api'
import { ErrorSnacks, ErrorBoundary } from './errors'
import searchQuantities from '../searchQuantities'
import '@testing-library/jest-dom/extend-expect' // Adds convenient expect-methods
import { keycloakBase } from '../config'
const fs = require('fs')
const crypto = require('crypto')
const { execSync } = require('child_process')
const keycloakURL = `${keycloakBase}realms/fairdi_nomad_test/protocol/openid-connect/token`

/*****************************************************************************/
// Renders

// Mocked keycloak setup. The real keycloak does not work within the test
// environment because access to the test realm is limited. Inspired by:
// https://stackoverflow.com/questions/63627652/testing-pages-secured-by-react-keycloak
let mockKeycloak = {
  init: jest.fn().mockResolvedValue(true),
  updateToken: jest.fn(),
  login: jest.fn(),
  logout: jest.fn(),
  register: jest.fn(),
  accountManagement: jest.fn(),
  createLoginUrl: jest.fn(),
  loadUserInfo: jest.fn(),
  authenticated: false,
  token: '',
  refreshToken: ''
}

const getRefreshToken = (username, password) => {
  let command = `curl -X POST ${keycloakURL} \\
    -H 'cache-control: no-cache' -H 'content-type: application/x-www-form-urlencoded' \\
    -d 'username=${username}&grant_type=password&password=${password}&client_id=nomad_gui_dev'`
  let response = require('child_process').execSync(command).toString()
  response = JSON.parse(response)
  if (response.error !== undefined) throw Error(response.error)
  mockKeycloak.login = (username ? () => login(username, password) : jest.fn())
  mockKeycloak.loadUserInfo = (username ? () => loadUserInfo(username) : jest.fn())
  return response
}

const updateToken = (refresh_token) => {
  return new Promise(async (resolve, reject) => {
    try {
      let command = `curl -X POST ${keycloakURL} \\
    -H 'cache-control: no-cache' -H 'content-type: application/x-www-form-urlencoded' \\
    -d 'refresh_token=${refresh_token}&grant_type=refresh_token&client_id=nomad_gui_dev'`
      let response = require('child_process').execSync(command).toString()
      response = JSON.parse(response)
      if (response.error !== undefined) return {}
      let authenticated = response.access_token !== undefined
      mockKeycloak.updateToken = (authenticated ? () => updateToken(response.refresh_token) : jest.fn())
      mockKeycloak.authenticated = authenticated
      mockKeycloak.token = (authenticated ? response.access_token : '')
      mockKeycloak.refreshToken = (authenticated ? response.refresh_token : '')
      resolve(response)
    } catch (error) {
      reject(new Error(error))
    }
  })
}

const testUsers = {
  'test': {
    'sub': '68878af7-6845-46c0-b2c1-250d4d8eb470',
    'email_verified': true,
    'name': 'Markus Scheidgen',
    'preferred_username': 'test',
    'given_name': 'Markus',
    'family_name': 'Scheidgen',
    'email': 'markus.scheidgen@fhi-berlin.de'
  },
  'scooper': {
    'sub': 'a03af8b6-3aa7-428a-b3b1-4a6317e576b6',
    'email_verified': true,
    'name': 'Sheldon Cooper',
    'preferred_username': 'scooper',
    'given_name': 'Sheldon',
    'family_name': 'Cooper',
    'email': 'sheldon.cooper@nomad-coe.eu'
  },
  'ttester': {
    'sub': '54cb1f64-f84e-4815-9ade-440ce0b5430f',
    'email_verified': true,
    'name': 'Test Tester',
    'preferred_username': 'ttester',
    'given_name': 'Test',
    'family_name': 'Tester',
    'email': 'test@nomad-coe.eu'
  },
  'admin': {
    'sub': 'c97facc2-92ec-4fa6-80cf-a08ed957255b',
    'email_verified': true,
    'name': 'Admin Administrator',
    'preferred_username': 'admin',
    'given_name': 'Admin',
    'family_name': 'Administrator',
    'email': 'markus.scheidgen@physik.hu-berlin.de'
  }
}

const loadUserInfo = (username) => {
  return new Promise(async (resolve, reject) => {
    try {
      const user = testUsers[username]
      resolve(user)
    } catch (error) {
      reject(new Error(error))
    }
  })
}

const login = (username, password) => {
  if ((username === undefined || username === '') && (password === undefined || password === '')) return
  let response = getRefreshToken(username, password)
  let authenticated = response.access_token !== undefined
  if (authenticated) updateToken(response.refresh_token)
}

const logout = () => {
  mockKeycloak.updateToken = jest.fn()
  mockKeycloak.authenticated = false
  mockKeycloak.token = ''
  mockKeycloak.refreshToken = ''
}

export const KeycloakProviderMock = (props) => {
  const { children } = props
  return (
    <KeycloakProvider keycloak={mockKeycloak}>
      {children}
    </KeycloakProvider>
  )
}

KeycloakProviderMock.propTypes = {
  children: PropTypes.node
}

const mockInitialized = true
jest.mock('react-keycloak', () => {
  const originalModule = jest.requireActual('react-keycloak')
  return {
    ...originalModule,
    useKeycloak: () => [
      mockKeycloak,
      mockInitialized
    ]
  }
})

/**
 * Default render function.
 */
export const WrapperDefault = ({children}) => {
  return <KeycloakProviderMock>
    <RecoilRoot>
      <APIProvider>
        <Router history={createBrowserHistory({basename: process.env.PUBLIC_URL})}>
          <MemoryRouter>
            <ErrorSnacks>
              <ErrorBoundary>
                {children}
              </ErrorBoundary>
            </ErrorSnacks>
          </MemoryRouter>
        </Router>
      </APIProvider>
    </RecoilRoot>
  </KeycloakProviderMock>
}

WrapperDefault.propTypes = {
  children: PropTypes.node
}

const renderDefault = (ui, options) =>
  render(ui, {wrapper: WrapperDefault, ...options})

export { renderDefault as render }

/**
 * Render function without API. Use for tests that do not require the API to
 * avoid error messages due to finishing tests before API info is retrieved.
 */
export const WrapperNoAPI = ({children}) => {
  return <KeycloakProviderMock>
    <RecoilRoot>
      <Router history={createBrowserHistory({basename: process.env.PUBLIC_URL})}>
        <MemoryRouter>
          <ErrorSnacks>
            <ErrorBoundary>
              {children}
            </ErrorBoundary>
          </ErrorSnacks>
        </MemoryRouter>
      </Router>
    </RecoilRoot>
  </KeycloakProviderMock>
}

WrapperNoAPI.propTypes = {
  children: PropTypes.node
}

export const renderNoAPI = (ui, options) =>
  render(ui, {wrapper: WrapperNoAPI, ...options})

/*****************************************************************************/
// Queries
/**
 * Query for finding an MUI MenuItem with the given text content.
 */
const queryAllByMenuItem = (c, value) =>
  queryAllByText(c, value, {selector: 'li'})
const multipleErrorMenuItem = (c, value) =>
  `Found multiple elements with the option value of: ${value}`
const missingErrorMenuItem = (c, value) =>
  `Unable to find an element with the option value of: ${value}`
const [
  queryByMenuItem,
  getAllByMenuItem,
  getByMenuItem,
  findAllByMenuItem,
  findByMenuItem
] = buildQueries(queryAllByMenuItem, multipleErrorMenuItem, missingErrorMenuItem)
const byMenuItem = {
  queryByMenuItem,
  queryAllByMenuItem,
  getByMenuItem,
  getAllByMenuItem,
  findAllByMenuItem,
  findByMenuItem
}

// Override default screen method by adding custom queries into it
const boundQueries = Object.entries(byMenuItem).reduce(
  (queries, [queryName, queryFn]) => {
    queries[queryName] = queryFn.bind(null, document.body)
    return queries
  },
  {}
)
const customScreen = { ...screen, ...boundQueries }
export { customScreen as screen }

/*****************************************************************************/
// Expects

/**
 * Used to find and test data displayed using the Quantity-component.
 *
 * @param {string} name Full metainfo name for the quantity.
 * @param {object} data A data object or a data value. If a data object is provided, the
 * metainfo name will be used to query for the final value.
 * @param {string} label Label to search for, optional, by default read using metainfo
 * name.
 * @param {string} description Description to search for, optional, by default read using
 * metainfo name.
 * @param {object} root The container to work on.
*/
export function expectQuantity(name, data, label = undefined, description = undefined, root = screen) {
  description = description || searchQuantities[name].description.replace(/\n/g, ' ')
  label = label || searchQuantities[name].name.replace(/_/g, ' ')
  const value = isPlainObject(data) ? get(data, name) : data
  const element = root.getByTitle(description)
  expect(root.getByText(label)).toBeInTheDocument()
  expect(within(element).getByText(value)).toBeInTheDocument()
}

/*****************************************************************************/
// Misc

/**
 * Used to prepare an API state for a test.
 *
 * Primarily uses a pre-recorded API snapshot file if one is available.
 * Otherwise records the API traffic into a file.
 *
 * @param {string} state Name of the state to prepare.
 * @param {string} path Path of the file in which the API traffic will be recorder or
 * from which it will be read.
 * @param {string} mode Snapshot usage mode:
 *   -'r': Only read
 *   -'w': Only write
 *   -'rw': Read and write
 */
let filepath
let responseCapture = {}
const readMode = process.env.READ_MODE || 'snapshot'
const writeMode = process.env.WRITE_MODE || 'none'
const configPath = 'nomad-test.yaml'
if (!fs.existsSync(`../${configPath}`)) {
  throw Error(`
    Could not find the NOMAD config file for testing at ../${configPath}. Note
    that the test environment should use a configuration that does not interfere
    with other NOMAD installations.
  `)
}
export function startAPI(state, path, username = '', password = '') {
  const jsonPath = `${path}.json`
  responseCapture[jsonPath] = {}

  // Prepare API state for reading responses directly from it.
  login(username, password)
  if (readMode === 'api') {
    // Prepare the test state
    const splits = state.split('.')
    const module = splits.slice(0, splits.length - 1).join('.')
    const func = splits[splits.length - 1]
    execSync(`
cd ..;
export NOMAD_CONFIG=${configPath};
python -c "
from ${module} import ${func}
${func}()"`)
  // Use API responses from snapshot files
  } else if (readMode === 'snapshot') {
    if (fs.existsSync(jsonPath)) {
      const apiSnapshot = JSON.parse(fs.readFileSync(jsonPath))
      const read = () => {
        const counterMap = {}
        return async (req, res, ctx) => {
          const hash = hashRequest(req)
          const counter = counterMap[hash] || 0
          const responseSnap = apiSnapshot[hash][counter].response
          const isJson = !responseSnap.headers['content-type'] || responseSnap.headers['content-type'] === 'application/json'
          const response = res(
            ctx.status(responseSnap.status),
            ctx.set(responseSnap.headers),
            isJson ? ctx.json(responseSnap.body) : ctx.body(responseSnap.body)
          )
          counterMap[hash] = counter + 1
          return response
        }
      }
      const readHandlers = [
        rest.get('*', read()),
        rest.post('*', read()),
        rest.patch('*', read()),
        rest.put('*', read()),
        rest.delete('*', read())
      ]
      server.use(...readHandlers)
    } else {
      throw Error(`Could not find the snapshot file: ${jsonPath}`)
    }
  }
  // Save results into a snapshot file using MSW. Notice that certain highly
  // static information will not be cpature, e.g. info.
  if (writeMode === 'snapshot') {
    filepath = jsonPath
    const capture = async (req, res, ctx) => {
      const hash = hashRequest(req)
      const response = await ctx.fetch(req)
      const headers = Object.fromEntries(response.headers.entries())
      delete headers['date']
      const isJson = !headers['content-type'] || headers['content-type'] === 'application/json'
      const body = isJson ? await response.json() : await response.text()
      const snapshot = {
        request: {
          url: req.url,
          method: req.method,
          body: req.body,
          headers: Object.fromEntries(req.headers.entries())
        },
        response: {
          status: response.status,
          body: body,
          headers: headers
        }
      }
      if (hash in responseCapture[filepath]) {
        responseCapture[filepath][hash].push(snapshot)
      } else {
        responseCapture[filepath][hash] = [snapshot]
      }
      return res(
        ctx.status(response.status),
        ctx.set(headers),
        isJson ? ctx.json(body) : ctx.body(body)
      )
    }
    const captureHandlers = [
      rest.get('*', capture),
      rest.post('*', capture),
      rest.patch('*', capture),
      rest.put('*', capture),
      rest.delete('*', capture)
    ]

    server.use(...captureHandlers)
  }
}

/**
 * Creates a hash for an HTTP request.
 */
function hashRequest(req) {
  const url = req.params['0']
  const method = req.method
  const body = req.body
  return crypto
    .createHash('md5')
    .update(url, method, body)
    .digest('hex')
}

/**
 * Used to close the API state.
 *
 * If an API response snapshot has been recorded, the snaphost will be saved as
 * a JSON file and the response handlers used to capture the API calls will be
 * cleared.
 */
export function closeAPI() {
  logout()

  // Tear down the test state when running a live API
  if (readMode === 'api') {
    execSync(`
cd ..;
export NOMAD_CONFIG=${configPath};
python -c "
from nomad import infrastructure
infrastructure.setup_mongo()
infrastructure.setup_elastic()
infrastructure.reset(True)"
  `)
  }
  // Write snapshot file
  if (writeMode === 'snapshot' && filepath) {
    fs.writeFile(filepath, JSON.stringify(responseCapture[filepath], null, 2), 'utf8', function(err) {
      if (err) {
        console.log('An error occured while writing JSON API response to file.')
        return console.log(err)
      }
    })
  }
  server.resetHandlers()
  responseCapture[filepath] = {}
  filepath = undefined
}

/**
 * Utility function for emulating delayed execution.
 *
 * @param {*} value value to return after delay
 * @param {number} ms delay in milliseconds
 */
export function wait(value, ms = 100) {
  return new Promise(resolve => setTimeout(() => resolve(value), ms))
}

/**
 * Utility function for getting the archive and search index from a JSON archive
 * file.
 *
 * @param {string} path Path to a JSON archive file.
 */
export async function readArchive(path) {
  const archive = await import(path)
  return [archive, {...archive, ...archive.metadata}]
}

// Map from entry_id to an archive
export const archives = new Map()
const archive = getArchive()
archives.set(archive.metadata.entry_id, archive)

/**
 * Convenience utility providing an easy way to find a button associated with a certain text,
 * without needing to care about the exact implementation (i.e. exactly how the text and the
 * button are related). We try a number of different ways to find a matching button. Note,
 * that if the button/text association is created in some unusual way, this method may not
 * find it. The role must be 'button' for this method to find it. If no button can be found,
 * null is returned.
 * @param {string} text Text to look for. Can be a string (case sensitive) or regex.
 * @param {*} findWithin An optional object within which to search for the button.
 */
export function findButton(text, findWithin = undefined) {
  const buttons = findWithin ? within(findWithin).queryAllByRole('button') : screen.queryAllByRole('button')
  for (const button of buttons) {
    if (typeof text === 'string') {
      // String provided - exact match required.
      if (button.title === text || button.text === text) {
        return button
      }
    } else {
      // Regex provided
      if (text.test(button.title) || text.test(button.text)) {
        return button
      }
    }
    if (within(button).queryByTitle(text) || within(button).queryByText(text)) {
      return button
    }
  }
  return null
}
