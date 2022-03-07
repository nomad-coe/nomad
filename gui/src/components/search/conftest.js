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
import PropTypes from 'prop-types'
import { screen, WrapperDefault, WrapperNoAPI } from '../conftest'
import { render } from '@testing-library/react'
import { SearchContext } from './SearchContext'
import { filterData } from './FilterRegistry'

/*****************************************************************************/
// Renders
/**
 * Entry search render.
 */
const WrapperSearch = ({children}) => {
  return <WrapperDefault>
    <SearchContext resource="entries">
      {children}
    </SearchContext>
  </WrapperDefault>
}

WrapperSearch.propTypes = {
  children: PropTypes.node
}

export const renderSearchEntry = (ui, options) =>
  render(ui, {wrapper: WrapperSearch, ...options})

/**
 * Entry search render without API.
 */
const WrapperNoAPISearch = ({children}) => {
  return <WrapperNoAPI>
    <SearchContext resource="entries">
      {children}
    </SearchContext>
  </WrapperNoAPI>
}

WrapperNoAPISearch.propTypes = {
  children: PropTypes.node
}

export const renderNoAPISearchEntry = (ui, options) =>
  render(ui, {wrapper: WrapperNoAPISearch, ...options})

/*****************************************************************************/
// Expects
/**
 * Tests that the initial state of an InputHeader is correct.
 *
 * @param {string} Quantity Full metainfo name for the quantity.
 * @param {boolean} disableScale Is the statistics scaling is disabled.
 * @param {object} root The container to work on.
 */
export function expectInputHeader(quantity, disableScale, root = screen) {
  const data = filterData[quantity]
  const label = data.label
  const description = data.description
  expect(root.getByText(label)).toBeInTheDocument()
  expect(root.getByTitle(description)).toBeInTheDocument()
  if (!disableScale) {
    const scale = data.scale
    expect(root.getByText(scale)).toBeInTheDocument()
  }
}
