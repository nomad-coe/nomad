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
import { get, isNil, isPlainObject } from 'lodash'
import { rest } from 'msw'
import PropTypes from 'prop-types'
import { RecoilRoot } from 'recoil'
import DateFnsUtils from '@date-io/date-fns'
import { MuiPickersUtilsProvider } from '@material-ui/pickers'
import {
  act,
  render,
  screen,
  within as originalWithin,
  buildQueries,
  queryAllByText,
  queryAllByRole,
  queries
} from '@testing-library/react'
import { prettyDOM, getDefaultNormalizer } from '@testing-library/dom'
import { seconds, server } from '../setupTests'
import { Router, MemoryRouter } from 'react-router-dom'
import { createBrowserHistory } from 'history'
import { APIProvider } from './api'
import { ErrorSnacks, ErrorBoundary } from './errors'
import DataStore from './DataStore'
import searchQuantities from '../searchQuantities'
import { filterData } from './search/FilterRegistry'
import { keycloakBase } from '../config'
import { useKeycloak } from '@react-keycloak/web'
import { GlobalMetainfo } from './archive/metainfo'

beforeEach(async () => {
  // For some strange reason, the useKeycloak mock gets reset if we set it earlier
  if (!useKeycloak()?.keycloak) {
    useKeycloak.mockReturnValue(
      {
        keycloak: {
          // Default Keycloak mock
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
        },
        initialized: true
      })
  }

  // Mock the window.ResizeObserver (jest does not seem to mock it in a way that works)
  class ResizeObserver {
    observe() {}
    unobserve() {}
    disconnect() {}
  }
  window.ResizeObserver = ResizeObserver
})

const fs = require('fs')
const crypto = require('crypto')
const { execSync } = require('child_process')
const keycloakURL = `${keycloakBase}realms/fairdi_nomad_test/protocol/openid-connect/token`

jest.mock('@react-keycloak/web', () => {
  const originalModule = jest.requireActual('@react-keycloak/web')

  return {
    __esModule: true,
    ...originalModule,
    useKeycloak: jest.fn()
  }
})

/*****************************************************************************/
// Renders
/**
 * Default render function.
 */
export const WrapperDefault = ({children}) => {
  return <RecoilRoot>
    <APIProvider>
      <MuiPickersUtilsProvider utils={DateFnsUtils}>
        <ErrorSnacks>
          <ErrorBoundary>
            <DataStore>
              <GlobalMetainfo>
                <Router history={createBrowserHistory({basename: process.env.PUBLIC_URL})}>
                  <MemoryRouter>
                    {children}
                  </MemoryRouter>
                </Router>
              </GlobalMetainfo>
            </DataStore>
          </ErrorBoundary>
        </ErrorSnacks>
      </MuiPickersUtilsProvider>
    </APIProvider>
  </RecoilRoot>
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
  return <RecoilRoot>
    <MuiPickersUtilsProvider utils={DateFnsUtils}>
      <Router history={createBrowserHistory({basename: process.env.PUBLIC_URL})}>
        <MemoryRouter>
          <ErrorSnacks>
            <ErrorBoundary>
              {children}
            </ErrorBoundary>
          </ErrorSnacks>
        </MemoryRouter>
      </Router>
    </MuiPickersUtilsProvider>
  </RecoilRoot>
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
  queryAllByMenuItem,
  queryByMenuItem,
  getAllByMenuItem,
  getByMenuItem,
  findAllByMenuItem,
  findByMenuItem
}

/**
 * Query for finding a DOM element with a certain role that is "associated" with a certain text.
 * The text must either be found as title or text ...
 *  - *on* the returned element itself
 *  - *within* the returned element itself, or
 *  - within an *ancestor* of the returned element. This requires that we can find an ancestor
 *    containing both A) the desired text and B) one and only one element which has the desired role,
 *    and C) that the ancestor is no more than ancestorDepth ancestors removed from the returned element.
 *    Setting ancestorDepth to 0 means no ancestors are considered, i.e. the text must be found
 *    on or within an element with the desired role.
 */
const queryAllByRoleAndText = (c, role, text, ancestorDepth = 5) => {
  const elementsWithRole = c.queryAllByRole ? c.queryAllByRole(role) : queryAllByRole(c, role)
  const rv = []
  for (const elementWithRole of elementsWithRole) {
    let element = elementWithRole
    for (let i = 0; i <= ancestorDepth; i++) {
      if (element !== elementWithRole && within(element).queryAllByRole(role).length > 1) {
        break // There must be only one element with the desired role within the node, or we stop.
      }
      if (typeof text === 'string') {
        // String provided
        if (element.title === text || element.text === text) {
          rv.push(elementWithRole)
          break
        }
      } else {
        // Regex provided
        if (text.test(element.title) || text.test(element.text)) {
          rv.push(elementWithRole)
          break
        }
      }
      if (within(element).queryByTitle(text) || within(element).queryByText(text)) {
        rv.push(elementWithRole)
        break
      }
      element = element.parentElement
      if (!element) {
        break
      }
    }
  }
  return rv
}
const multipleErrorRoleAndText = (c, role, text) =>
  `Found multiple elements with role '${role}' and text '${text}'`
const missingErrorRoleAndText = (c, role, text) =>
  `Unable to find an element with role '${role}' and text '${text}'`
const [
  queryByRoleAndText,
  getAllByRoleAndText,
  getByRoleAndText,
  findAllByRoleAndText,
  findByRoleAndText
] = buildQueries(queryAllByRoleAndText, multipleErrorRoleAndText, missingErrorRoleAndText)
const byRoleAndText = {
  queryAllByRoleAndText,
  queryByRoleAndText,
  getAllByRoleAndText,
  getByRoleAndText,
  findAllByRoleAndText,
  findByRoleAndText
}

/**
 * Query for finding buttons associated with a given text somehow. The text may be a string
 * or a regex. An element will match if it has role = "button" and, for example, the text
 * matches the button's title or is found somewhere within the button. The goal is for this
 * method to handle all types of buttons we use, and it should always return an object which
 * is of a button type and therefore can be used for firing events etc.
 */
const queryAllByButtonText = (c, text, ancestorDepth = 0) => {
  return queryAllByRoleAndText(c, 'button', text, ancestorDepth)
}
const multipleErrorButtonText = (c, text) =>
  `Found multiple buttons with the text: ${text}`
const missingErrorButtonText = (c, text) =>
  `Unable to find a button with the text: ${text}`
const [
  queryByButtonText,
  getAllByButtonText,
  getByButtonText,
  findAllByButtonText,
  findByButtonText
] = buildQueries(queryAllByButtonText, multipleErrorButtonText, missingErrorButtonText)
const byButtonText = {
  queryAllByButtonText,
  queryByButtonText,
  getAllByButtonText,
  getByButtonText,
  findAllByButtonText,
  findByButtonText
}

// Override default screen method by adding custom queries into it
const customQueries = {...byMenuItem, ...byRoleAndText, ...byButtonText}
const defaultAndCustomQueries = {...queries, ...customQueries}

const boundCustomQueries = Object.entries(customQueries).reduce(
  (queries, [queryName, queryFn]) => {
    queries[queryName] = queryFn.bind(null, document.body)
    return queries
  },
  {}
)
const customScreen = { ...screen, ...boundCustomQueries }
export { customScreen as screen }

export function within(element, queriesToBind = defaultAndCustomQueries) {
  return originalWithin(element, queriesToBind)
}

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
  description = description || searchQuantities[name].description
  label = label || filterData?.[name]?.label?.toLowerCase() || searchQuantities[name].name.replace(/_/g, ' ')
  const value = isPlainObject(data) ? get(data, name) : data
  const element = root.getByTitle(description, {normalizer: getDefaultNormalizer({trim: false, collapseWhitespace: false})})
  expect(root.getByText(label)).toBeInTheDocument()
  expect(within(element).getByText(value)).toBeInTheDocument()
}

/*****************************************************************************/
// Misc
let filepath
const responseCapture = {}
const readMode = process.env.READ_MODE || 'snapshot'
const writeMode = process.env.WRITE_MODE || 'none'
const configPath = 'gui/tests/nomad.yaml'
if (!fs.existsSync(`../${configPath}`)) {
  throw Error(`
    Could not find the NOMAD config file for testing at ../${configPath}. Note
    that the test environment should use a configuration that does not interfere
    with other NOMAD installations.
  `)
}
/**
 * Used to prepare an API state for a test.
 *
 * Primarily uses a pre-recorded API snapshot file if one is available.
 * Otherwise records the API traffic into a file.
 *
 * @param {string} state Name of the state to prepare.
 * @param {string} path Path of the file in which the API traffic will be recorded or
 * from which it will be read.
 * @param {string} username Username to login with.
 * @param {string} password Password for the username.
 */
export async function startAPI(state, path, username = '', password = '') {
  await mockKeycloak(username, password)

  // Prepare API state for reading responses directly from it.
  const jsonPath = `${path}.json`
  responseCapture[jsonPath] = {}
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
    // We need to wait here for some tasks to finish. TODO: resolve why this is
    // needed, and how could this simple wait be replaced with something
    // better.
    await wait(undefined, 3 * seconds)
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
          url: req.url.toString(),
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
 * Creates a mocked instance of useKeycloak for the given user. The real
 * keycloak does not work within the test environment because access to the test
 * realm is limited. Inspired by:
 * https://stackoverflow.com/questions/63627652/testing-pages-secured-by-react-keycloak
 */
async function mockKeycloak(username, password) {
  const login = async (username, password) => {
    if ((username === undefined || username === '') && (password === undefined || password === '')) return
    const response = getRefreshToken(username, password)
    const authenticated = response.access_token !== undefined
    if (authenticated) await updateToken(response.refresh_token)
  }

  const logout = () => {
    mockedKeycloak.updateToken = jest.fn()
    mockedKeycloak.authenticated = false
    mockedKeycloak.token = ''
    mockedKeycloak.refreshToken = ''
  }

  const getRefreshToken = (username, password) => {
    const command = `curl -s -X POST ${keycloakURL} \\
      -H 'cache-control: no-cache' -H 'content-type: application/x-www-form-urlencoded' \\
      -d 'username=${username}&grant_type=password&password=${password}&client_id=nomad_gui_dev'`
    let response = require('child_process').execSync(command).toString()
    response = JSON.parse(response)
    if (response.error !== undefined) throw Error(response.error)
    mockedKeycloak.login = (username ? () => login(username, password) : jest.fn())
    mockedKeycloak.loadUserInfo = (username ? () => loadUserInfo(username) : jest.fn())
    return response
  }

  const updateToken = (refresh_token) => {
    return new Promise((resolve, reject) => {
      try {
        const command = `curl -s -X POST ${keycloakURL} \\
      -H 'cache-control: no-cache' -H 'content-type: application/x-www-form-urlencoded' \\
      -d 'refresh_token=${refresh_token}&grant_type=refresh_token&client_id=nomad_gui_dev'`
        let response = require('child_process').execSync(command).toString()
        response = JSON.parse(response)
        if (response.error !== undefined) return {}
        const authenticated = response.access_token !== undefined
        mockedKeycloak.updateToken = (authenticated ? () => updateToken(response.refresh_token) : jest.fn())
        mockedKeycloak.authenticated = authenticated
        mockedKeycloak.token = (authenticated ? response.access_token : '')
        mockedKeycloak.refreshToken = (authenticated ? response.refresh_token : '')
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
    return new Promise((resolve, reject) => {
      try {
        const user = testUsers[username]
        resolve(user)
      } catch (error) {
        reject(new Error(error))
      }
    })
  }

  const mockedKeycloak = {
    init: jest.fn().mockResolvedValue(true),
    updateToken: updateToken,
    login: login,
    logout: logout,
    register: jest.fn(),
    accountManagement: jest.fn(),
    createLoginUrl: jest.fn(),
    loadUserInfo: loadUserInfo,
    authenticated: false,
    token: '',
    refreshToken: ''
  }

  if (username && password) {
    await login(username, password)
  }

  useKeycloak.mockReturnValue({keycloak: mockedKeycloak, initialized: true})
}

/**
 * Creates a hash for an HTTP request.
 */
function hashRequest(req) {
  const url = req.url.toString()
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
  // Tear down the test state when running a live API
  if (readMode === 'api') {
    execSync(`
cd ..;
export NOMAD_CONFIG=${configPath};
python -c "
import time
from nomad import infrastructure
infrastructure.setup_mongo()
infrastructure.setup_elastic()
infrastructure.reset(True)"`)
  }
  // Write snapshot file
  if (writeMode === 'snapshot') {
    if (isNil(filepath)) {
      throw Error(
        'The snapshot filepath was not specified. Did you remember to use "await' +
        ' startAPI" to wait for the startup to finish correctly?'
      )
    }
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

  // Restore the default keycloak mock that has no user credentials
  useKeycloak.mockRestore()
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
 * Utility function for waiting for the GUI within an "act" statement (so that jest does
 * not print warnings if there are state updates while waiting). Intended to be used as
 * a LAST RESORT when updates are expected in the GUI and there are either no good ways of
 * determining when they are completed, or we use this as a temporary workaround until we
 * have a proper solution.
 * @param {number} ms delay in milliseconds
 * @param waitInActualTest force to wait in the actual test e.g. when the delay comes from a debounce
 */
export async function waitForGUI(ms = 1000, waitInActualTest = false) {
  if (process.env.WAIT_FOR_GUI !== 'None' || waitInActualTest) {
    await act(async () => { await new Promise(resolve => setTimeout(resolve, ms)) })
  }
}

/**
 * Utility function for getting the archive and search index from a JSON archive
 * file.
 *
 * @param {string} path Path to a JSON archive file.
 */
export async function readArchive(path) {
  const archive = await import(path)
  const index = {...archive, ...archive.metadata}
  const properties = new Set(index?.results
    ? index.results.properties.available_properties
    : [])
  return {archive, index, properties}
}

export const consoleSpies = {}
const consoleIgnoreStrings = [
  'Warning: You seem to have overlapping act()',
  'Warning: An update to %s inside a test was not wrapped in act'
]

/**
 * Utility for spying on the console, and yielding errors if something is printed to it
 * (since most of the time, we don't want anything to be printed to the console).
 * Typical usage: call blockConsoleOutput before the test, and unblockConsoleOutput after,
 * for example using beforeEach and afterEach. Note, some common jest warnings are ignored.
 */
export function blockConsoleOutput() {
  if (consoleSpies.logSpy) {
    throw Error('Already spying on the console output')
  }
  consoleSpies.logSpy = jest.spyOn(console, 'log')
  consoleSpies.errorSpy = jest.spyOn(console, 'error')
}

/**
 * Returns a list with the strings printed to the console, except excluded jest warnings.
 */
export function filteredConsoleOutput() {
  const rv = []
  for (const consoleSpy of [consoleSpies.logSpy, consoleSpies.errorSpy]) {
    for (const call of consoleSpy.mock.calls) {
      const message = '' + call[0]
      let isOk = false
      for (const s of consoleIgnoreStrings) {
        if (message.startsWith(s)) {
          isOk = true
          break
        }
      }
      if (!isOk) {
        rv.push(message)
      }
    }
  }
  return rv
}

/**
 * Expects that nothing has been written to the console (so far). Note, some common jest warnings are ignored.
 */
export function expectNoConsoleOutput() {
  if (!consoleSpies.logSpy) {
    throw Error('Need to call blockConsoleOutput before using this method!')
  }
  const consoleOutput = filteredConsoleOutput()
  if (consoleOutput.length) {
    for (const message of consoleOutput) {
      process.stdout.write(`Unexpected console output: ${message}\n`)
    }
    expect(consoleOutput.length).toBe(0) // Fails
  }
}

/**
 * Checks that nothing has been written to the console, and removes the log spies
 * (see blockConsoleOutput)
 */
export function unblockConsoleOutput() {
  try {
    expectNoConsoleOutput()
  } finally {
    consoleSpies.logSpy.mockRestore()
    consoleSpies.errorSpy.mockRestore()
    consoleSpies.logSpy = consoleSpies.errorSpy = null
  }
}

/**
 * Utility for dumping a testing-library dom object to file, for inspection.
 * @param {*} domObj
 * @param {*} filePath Path to file. Default is 'domdump.html'
 */
export function dumpPrettyDom(domObj, filePath = 'domdump.html') {
  return new Promise((resolve, reject) => {
    fs.writeFile(filePath, prettyDOM(domObj, 10000000, {highlight: false}), 'utf8', (err) => {
      if (err) reject(err)
      else resolve()
    })
  })
}
