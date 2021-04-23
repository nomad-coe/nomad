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
import { compose } from 'recompose'
import { render } from '@testing-library/react'
import { apiContext as apiContextV0 } from './components/api'
import { apiContext as apiContextV1 } from './components/apiv1'
import { Router } from 'react-router-dom'
import { archiveDftBulk } from '../tests/archive'
import history from './history'

/**
 * Mock for useKeycloak in the package 'react-keycloak'. Notice that this is
 * for Keycloak 9.x which has a slighly different API from the newer versions
 * (e.g. loadUserInfo which does not return promises)
 */
export function useKeycloakMock() {
  const token = 'A random string that is non zero length'
  const userProfile = { username: 'test', email: 'test@testdomain.com', firstName: 'Test', lastName: 'User' }
  const realmAccess = { roles: ['admin', 'auditor', 'user'] }
  let authenticated = false

  const authClient = {
    authenticated,
    hasRealmRole(role) {
      return true
    },
    hasResourceRole(role) {
      return true
    },
    idToken: token,
    initialized: true,
    loadUserInfo() {
      return {success: (fn) => {
        fn(userProfile)
        return {error: (fn) => {}}
      }}
    },
    updateToken() {
      return {success: (fn) => {
        fn()
        return {error: (fn) => {}}
      }}
    },
    login() {
      authenticated = true
    },
    logout() {
      authenticated = false
    },
    profile: userProfile,
    realm: 'test_realm',
    realmAccess,
    refreshToken: token,
    token
  }
  return { initialized: true, keycloak: authClient }
}

/**
 * HOC for injecting a mocked Keycloak setup. Mocks the 'withKeycloak' HOC in
 * the package 'react-keycloak'
 */
export function withKeycloakMock(Component) {
  return function WrappedComponent(props) {
    const { keycloak, initialized } = useKeycloakMock()
    return <Component {...props} keycloakInitialized={initialized} keycloak={keycloak} />
  }
}

/**
 * HOC for injecting a mocked API implementation.
 */
export function withApiV0Mock(Component) {
  const apiValue = {
    api: {
      archive: (upload_id, calc_id) => {
        return wait(archiveDftBulk)
      }
    }
  }
  return <apiContextV0.Provider value={apiValue}>
    {Component}
  </apiContextV0.Provider>
}

/**
 * HOC for injecting a mocked API implementation.
 */
export function withApiV1Mock(Component) {
  const apiValue = {
    api: {
      results: (entry_id) => {
        return wait(archiveDftBulk)
      }
    }
  }
  return <apiContextV1.Provider value={apiValue}>
    {Component}
  </apiContextV1.Provider>
}

/**
 * HOC for Router dependency injection
 */
export function withRouterMock(Component) {
  return <Router history={history}>
    {Component}
  </Router>
}

/**
 * Utility function for rendering components that require a Router and an API.
 */
export function renderWithAPIRouter(component) {
  return render(compose(
    withRouterMock,
    withApiV0Mock,
    withApiV1Mock
  )(component))
}

/**
 * Utility function for emulating delayed execution.
 *
 * @param {*} value value to return after delay
 * @param {number} ms delay in milliseconds
 */
function wait(value, ms = 200) {
  return new Promise(resolve => setTimeout(() => resolve(value), ms))
}
