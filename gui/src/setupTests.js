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

import { setupServer } from 'msw/node'
import { configure } from '@testing-library/react'
import "../tests/env"
import "../tests/artifacts"
import '@testing-library/jest-dom' // Adds convenient expect-methods

/**
 * Code to configure or set up the testing framework before each test file in
 * the suite is executed. Contains e.g. global setup/teardown functionality for
 * tests.
 */

export const seconds = 1000
export const minutes = 60 * seconds
// Increased the default jest timeout for individual tests
// eslint-disable-next-line no-undef
jest.setTimeout(2 * minutes)

// Changes the default timeout of waitFor, find*-queries etc. We set this very generously,
// because when tests are run in parallel, something which normally renders quickly may take
// long to render if it competes with other tests over limited CPU resources.
configure({ asyncUtilTimeout: 15 * seconds })

// Mocks required by Plotly.js:
// https://github.com/plotly/react-plotly.js/issues/115
// https://stackoverflow.com/questions/48828759/unit-test-raises-error-because-of-getcontext-is-not-implemented
window.URL.createObjectURL = function() {}
HTMLCanvasElement.prototype.getContext = () => {}

const { ResizeObserver } = window
export const server = setupServer()

beforeAll(() => {
  // Start MSW server. Unhandled requests are bypassed without warning.
  server.listen({onUnhandledRequest: 'bypass'})

  // ResizeObserver mock init
  delete window.ResizeObserver
  window.ResizeObserver = jest.fn().mockImplementation(() => ({
    observe: jest.fn(),
    unobserve: jest.fn(),
    disconnect: jest.fn()
  }))
})

// This is a default resize-observer result. You should overwrite this if you
// expect different size in the tests.
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {width: 500, height: 500, ref: undefined}
  }}
})

afterEach(() => {
})

afterAll(() => {
  // Close MSW server
  server.close()

  // ResizeObserver mock cleanup
  window.ResizeObserver = ResizeObserver
  jest.restoreAllMocks()
})
