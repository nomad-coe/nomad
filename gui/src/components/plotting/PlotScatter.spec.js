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
import React, { useRef } from 'react'
import { renderNoAPI, screen } from '../conftest.spec'
import PlotScatter from './PlotScatter'
import { DType } from '../../utils'

const ScatterPlotWrapper = React.memo((props) => {
  const canvas = useRef()
  return <PlotScatter {...props} ref={canvas}/>
})

test.each([
  [
    'linear timestamp',
    {x: ["2014-12-01T00:00:00+00:00"], y: ["2015-12-01T00:00:00+00:00"]},
    {title: 'Test', dtype: DType.Timestamp, scale: 'linear'},
    ['Dec 1, 2014', 'Dec 1, 2015']
  ],
  [
    'log float',
    {x: [1, 10], y: [100, 1000]},
    {title: 'Test', dtype: DType.Float, scale: 'log'},
    ['1', '10', '100', '1000']
  ]
])('%s', async (id, data, axis, ticks) => {
  renderNoAPI(<ScatterPlotWrapper
    data={data}
    xAxis={axis}
    yAxis={axis}
  />)

  // Check that the correct plot ticks are found
  for (const tick of ticks) {
    expect(await screen.findByText(tick)).toBeInTheDocument()
  }
})
