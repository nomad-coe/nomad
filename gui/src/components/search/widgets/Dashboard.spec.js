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
import { render, screen, startAPI, closeAPI } from '../../conftest.spec'
import userEvent from '@testing-library/user-event'
import Dashboard from './Dashboard'
import { defaultFilterData } from '../FilterRegistry'
import { SearchContext } from '../SearchContext'
import {
  expectWidgetTerms,
  expectWidgetScatterPlot,
  expectInputHeader,
  expectInputRange
} from '../conftest.spec'

// Resize-detector size is adjusted for the widgets.
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {width: 1000, height: 300, ref: undefined}
  }}
})

describe('displaying an initial widget and removing it', () => {
  beforeAll(async () => {
    await startAPI('tests.states.search.search', 'tests/data/search/dashboard')
  })
  afterAll(() => closeAPI())

  test.each([
    [
      'terms',
      {
        type: 'terms',
        search_quantity: 'results.material.structural_type',
        scale: 'linear',
        editing: false,
        visible: true,
        layout: {
          sm: {x: Infinity, y: 0, w: 6, h: 9},
          md: {x: Infinity, y: 0, w: 6, h: 9},
          lg: {x: Infinity, y: 0, w: 6, h: 9},
          xl: {x: Infinity, y: 0, w: 6, h: 9},
          xxl: {x: Infinity, y: 0, w: 6, h: 9}
        }
      },
      async (widget, loaded) => await expectWidgetTerms(
        widget,
        loaded,
        ['2D', 'bulk', 'molecule / cluster'],
        'all'
      )
    ],
    [
      'histogram',
      {
        type: 'histogram',
        title: 'Test title',
        x: {search_quantity: 'results.material.n_elements'},
        y: {scale: 'linear'},
        editing: false,
        visible: true,
        layout: {
          sm: {x: Infinity, y: 0, w: 12, h: 9},
          md: {x: Infinity, y: 0, w: 12, h: 9},
          lg: {x: Infinity, y: 0, w: 12, h: 9},
          xl: {x: Infinity, y: 0, w: 12, h: 9},
          xxl: {x: Infinity, y: 0, w: 12, h: 9}
        }
      },
      async (widget, loaded) => await expectInputRange(widget, loaded, true, true)
    ],
    [
      'scatter_plot',
      {
        type: 'scatter_plot',
        title: 'Test title',
        x: {search_quantity: 'results.properties.optoelectronic.solar_cell.open_circuit_voltage'},
        y: {search_quantity: 'results.properties.optoelectronic.solar_cell.efficiency'},
        markers: {color: {search_quantity: 'results.properties.optoelectronic.solar_cell.short_circuit_current_density'}},
        size: 1000,
        autorange: true,
        editing: false,
        visible: true,
        layout: {
          sm: {x: Infinity, y: 0, w: 9, h: 6},
          md: {x: Infinity, y: 0, w: 9, h: 6},
          lg: {x: Infinity, y: 0, w: 9, h: 6},
          xl: {x: Infinity, y: 0, w: 9, h: 6},
          xxl: {x: Infinity, y: 0, w: 9, h: 6}
        }
      },
      async (widget, loaded) => await expectWidgetScatterPlot(widget, loaded)
    ],
    [
      'periodic_table',
      {
        type: 'periodic_table',
        search_quantity: 'results.material.elements',
        scale: 'linear',
        editing: false,
        visible: true,
        layout: {
          sm: {x: Infinity, y: 0, w: 12, h: 9},
          md: {x: Infinity, y: 0, w: 12, h: 9},
          lg: {x: Infinity, y: 0, w: 12, h: 9},
          xl: {x: Infinity, y: 0, w: 12, h: 9},
          xxl: {x: Infinity, y: 0, w: 12, h: 9}
        }
      },
      async (widget, loaded) => {
        await expectInputHeader(widget.search_quantity)
        // TODO: For some reason this test works out fine locally, but always
        // fails in the CI. Could not be resolved even by making the tests wait
        // much longer.
        // await expectPeriodicTable(
        //   widget.quantity,
        //   loaded,
        //   ['H', 'C', 'N', 'I', 'Pb', 'Ti', 'Zr', 'Nb', 'Hf', 'Ta']
        // )
      }
    ]
  ])('%s', async (type, widget, test) => {
    const initialDashboard = {widgets: [widget]}
    render(
      <SearchContext resource="entries" initialDashboard={initialDashboard}>
        <Dashboard />
      </SearchContext>
    )

    // Assert that the item is displayed initially
    await test(widget, false)

    // Remove widget, check that it is gone. A test id is used to fetch the
    // remove button since it is an icon that may appear in several locations.
    const removeButton = screen.getByTestId(`0-remove-widget`)
    const label = widget.title || defaultFilterData[widget.search_quantity].label
    expect(screen.queryByText(label, {exact: false})).toBeInTheDocument()
    await userEvent.click(removeButton)
    expect(screen.queryByText(label, {exact: false})).not.toBeInTheDocument()

    // TODO: Re-add item, check that it appears. This time the data is already
    // loaded.
  })
})
