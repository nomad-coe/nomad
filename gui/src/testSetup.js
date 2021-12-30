/* eslint-disable import/export */
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

/**
 * Code to configure or set up the testing environment. Will be executed in the
 * testing environment immediately before executing the test code itself. Notice
 * that these will execute before setupTests.js and here we can't e.g. define
 * test setup/teardown.
 */
import React from 'react'
import PropTypes from 'prop-types'
import { RecoilRoot } from 'recoil'
import { render } from '@testing-library/react'
import { Router, MemoryRouter } from 'react-router-dom'
import { getArchive } from '../tests/DFTBulk'
import { createBrowserHistory } from 'history'

// This extends the expect method with a lot of useful stuff, see:
// https://github.com/testing-library/jest-dom
import '@testing-library/jest-dom/extend-expect'

// Map from entry_id to an archive
export const archives = new Map()
const archive = getArchive()
archives.set(archive.metadata.entry_id, archive)

/**
 * Provides mocked App infrastructure for testing
 */
const AllTheProviders = ({children}) => {
  return <RecoilRoot>
    <Router history={createBrowserHistory({basename: process.env.PUBLIC_URL})}>
      <MemoryRouter>
        {children}
      </MemoryRouter>
    </Router>
  </RecoilRoot>
}
AllTheProviders.propTypes = {
  children: PropTypes.node
}

const customRender = (ui, options) =>
  render(ui, {wrapper: AllTheProviders, ...options})

// Re-export everything
export * from '@testing-library/react'

/**
 * Utility function for emulating delayed execution.
 *
 * @param {*} value value to return after delay
 * @param {number} ms delay in milliseconds
 */
export function wait(value, ms = 100) {
  return new Promise(resolve => setTimeout(() => resolve(value), ms))
}

// Override render method
export { customRender as render }
