/* eslint-disable */
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
import { get, isPlainObject } from 'lodash'
import { screen, within } from './testSetup'
import searchQuantities from './searchQuantities'

/**
 * Used to find and test data displayed using the Quantity-component.
 * 
 * @param name Full metainfo name for the quantity.
 * @param data A data object or a data value. If a data object is provided, the
 * metainfo name will be used to query for the final value.
 * @param label Label to search for, optional, by default read using metainfo
 * name.
 * @param label Description to search for, optional, by default read using
 * metainfo name.
*/
export function expectQuantity(name, data, label = undefined, description = undefined) {
  description = description || searchQuantities[name].description.replace(/\n/g, ' ')
  label = label || searchQuantities[name].name.replace(/_/g, ' ')
  const value = isPlainObject(data) ? get(data, name) : data
  const element = screen.getByTitle(description)
  expect(screen.getByText(label)).toBeInTheDocument()
  expect(within(element).getByText(value)).toBeInTheDocument()
}

export function expectPlotButtons(plot) {
  expect(within(plot).getByRole('button', {name: 'Reset view'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'Toggle fullscreen'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'Capture image'})).toBeInTheDocument()
  expect(within(plot).getByRole('button', {name: 'View data in the archive'})).toBeInTheDocument()
}

