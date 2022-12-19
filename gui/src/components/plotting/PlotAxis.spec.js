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
import { format, getTime } from 'date-fns'
import { render, screen } from '../conftest.spec'
import PlotAxis from './PlotAxis'
import { DType } from '../../utils'

const mockHeight = 50
const mockWidth = 300
const labelHeight = 12
const labelWidth = 50
const yearFormat = 'yyyy'
const monthFormat = 'MMM'
const dayFormat = 'MMM d'
const hourFormat = 'k:mm'
const minuteFormat = 'k:mm'
const secondFormat = 'k:mm:ss'

// Resize-detector size is adjusted to test the axis behaviour.
jest.mock('react-resize-detector', () => {
  return {useResizeDetector: () => {
    return {height: mockHeight, width: mockWidth, ref: undefined}
  }}
})

function getTs(...args) {
  return getTime(new Date(...args))
}

test.each([
  ['only two ticks', 'left', 'linear', 0, 1, DType.Int, 'SI', 2, ['0', '1']],
  ['do not show fractional ticks for integers', 'left', 'linear', 0, 1, DType.Int, 'SI', 6, ['0', '1']],
  ['show fractional ticks for floats', 'left', 'linear', 0, 1, DType.Float, 'scientific', 6, ['0', '0.5', '1']],
  ['do not not show more ticks than fit vertically', 'left', 'linear', 0, 100, DType.Int, 'SI', 20, ['0', '50', '100']],
  ['do not not show more ticks than fit horizontally', 'bottom', 'linear', 0, 100, DType.Int, 'SI', 20, ['0', '20', '40', '60', '80', '100']],
  ['thousands SI format', 'left', 'linear', 0, 5000, DType.Int, 'SI', 3, ['0', '2.5k', '5k']],
  ['thousands scientific format', 'left', 'linear', 0, 5000, DType.Int, 'scientific', 3, ['0', '2.5e+3', '5e+3']],
  ['millions SI format', 'left', 'linear', 0, 5000000, DType.Int, 'SI', 3, ['0', '2.5M', '5M']],
  ['millions scientific format', 'left', 'linear', 0, 5000000, DType.Int, 'scientific', 3, ['0', '2.5e+6', '5e+6']],
  ['negative numbers', 'left', 'linear', -5000, 0, DType.Int, 'SI', 3, ['-5k', '-2.5k', '0']],
  ['timestamp seconds', 'bottom', 'linear', 0, 2000, DType.Timestamp, undefined, 3, [0, 1000, 2000].map(x => format(x, secondFormat))],
  ['timestamp fifteen seconds', 'bottom', 'linear', 0, 30000, DType.Timestamp, undefined, 3, [0, 15000, 30000].map(x => format(x, secondFormat))],
  ['timestamp thirty seconds', 'bottom', 'linear', 0, 60000, DType.Timestamp, undefined, 3, [0, 30000, 60000].map(x => format(x, secondFormat))],
  ['timestamp minutes', 'bottom', 'linear', 0, 120000, DType.Timestamp, undefined, 3, [0, 60000, 120000].map(x => format(x, minuteFormat))],
  ['timestamp fifteen minutes', 'bottom', 'linear', 0, 1800000, DType.Timestamp, undefined, 3, [0, 900000, 1800000].map(x => format(x, minuteFormat))],
  ['timestamp thirty minutes', 'bottom', 'linear', 0, 3600000, DType.Timestamp, undefined, 3, [0, 1800000, 3600000].map(x => format(x, minuteFormat))],
  ['timestamp hours', 'bottom', 'linear', 0, 7200000, DType.Timestamp, undefined, 3, [0, 3600000, 7200000].map(x => format(x, hourFormat))],
  ['timestamp six hours', 'bottom', 'linear', getTs(1970, 0, 0, 3), getTs(1970, 0, 0, 15), DType.Timestamp, undefined, 2, ['6:00', '12:00']],
  ['timestamp twelve hours', 'bottom', 'linear', getTs(1970, 0, 0, 6), getTs(1970, 0, 1, 6), DType.Timestamp, undefined, 2, ['12:00', '24:00']],
  [
    'timestamp days',
    'bottom',
    'linear',
    getTs(1970, 0, 0, 12),
    getTs(1970, 0, 3, 12),
    DType.Timestamp,
    undefined,
    3,
    [getTs(1970, 0, 1), getTs(1970, 0, 2), getTs(1970, 0, 3)].map(x => format(x, dayFormat))],
  [
    'timestamp weeks',
    'bottom',
    'linear',
    getTs(1970, 0, 1, 12),
    getTs(1970, 0, 22, 12),
    DType.Timestamp,
    undefined,
    3,
    [getTs(1970, 0, 4, 12), getTs(1970, 0, 11, 12), getTs(1970, 0, 18, 12)].map(x => format(x, dayFormat))
  ],
  [
    'timestamp months',
    'bottom',
    'linear',
    getTs(1970, 0, 15, 12),
    getTs(1970, 3, 15, 12),
    DType.Timestamp,
    undefined,
    3,
    [getTs(1970, 1, 1, 12), getTs(1970, 2, 1, 12), getTs(1970, 3, 1, 12)].map(x => format(x, monthFormat))
  ],
  [
    'timestamp quarter years',
    'bottom',
    'linear',
    getTs(1970, 2, 1),
    getTs(1973, 10, 1),
    DType.Timestamp,
    undefined,
    3,
    [getTs(1971, 4, 1), getTs(1971, 7, 1), getTs(1972, 10, 1)].map(x => format(x, yearFormat))
  ],
  [
    'timestamp years',
    'bottom',
    'linear',
    getTs(1970, 6, 1),
    getTs(1973, 6, 1),
    DType.Timestamp,
    undefined,
    3,
    [getTs(1971, 0, 1, 12), getTs(1971, 0, 1, 12), getTs(1972, 0, 1, 12)].map(x => format(x, yearFormat))
  ]
])('%s', async (msg, placement, scale, min, max, dtype, mode, nLabels, labels) => {
  render(<PlotAxis
    min={min}
    max={max}
    labels={nLabels}
    labelWidth={placement === 'left' ? undefined : labelWidth}
    labelHeight={labelHeight}
    mode={mode}
    dtype={dtype}
    scale={scale}
    placement={placement}
  />)

  // Check that the correct labels are found
  for (const label of labels) {
    expect(await screen.findByText(label)).toBeInTheDocument()
  }

  // Check that no extra labels are found
  const allLabels = screen.queryAllByText(/.+/)
  expect(allLabels.length).toBe(labels.length)
})
