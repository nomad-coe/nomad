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

/**
 * Code to configure or set up the testing framework before each test file in
 * the suite is executed. Contains e.g. global setup/teardown functionality for
 * tests.
 */
global.nomadEnv = {
  'keycloakBase': 'https://nomad-lab.eu/fairdi/keycloak/auth/',
  // Use the production API
  // 'keycloakRealm': 'fairdi_nomad_prod',
  // 'keycloakClientId': 'nomad_public',
  // 'appBase': 'https://nomad-lab.eu/prod/v1',
  // Use the local API
  'keycloakRealm': 'fairdi_nomad_test',
  'keycloakClientId': 'nomad_gui_dev',
  'appBase': 'http://localhost:8000/fairdi/nomad/latest',
  'encyclopediaBase': 'https://nomad-lab.eu/prod/rae/encyclopedia/#',
  'debug': false,
  'version': {
    'label': '1.1.0',
    'isBeta': false,
    'isTest': true,
    'usesBetaData': true,
    'officialUrl': 'https://nomad-lab.eu/prod/rae/gui'
  },
  'aitoolkitEnabled': false,
  'oasis': false,
  'servicesUploadLimit': 10
}
// Increased the default jest timeout for individual tests
// eslint-disable-next-line no-undef
jest.setTimeout(60000)

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

afterEach(() => {
})

afterAll(() => {
  // Close MSW server
  server.close()

  // ResizeObserver mock cleanup
  window.ResizeObserver = ResizeObserver
  jest.restoreAllMocks()
})
