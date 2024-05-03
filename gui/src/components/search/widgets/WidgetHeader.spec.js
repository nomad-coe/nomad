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
import { renderSearchEntry, expectFilterTitle } from '../conftest.spec'
import WidgetHeader from './WidgetHeader'

test.each([
  [
    'label without unit',
    { quantity: 'results.material.n_elements' },
    { quantity: 'results.material.n_elements' }
  ],
  [
    'label with unit',
    { quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective'},
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      unit: 'eV'
    }
  ],
  [
    'full label',
    {quantity: 'results.properties.electronic.band_structure_electronic.band_gap.value'},
    {
      quantity: 'results.properties.electronic.band_structure_electronic.band_gap.value',
      label: 'Value',
      unit: 'eV'
    }
  ]
])('%s', async (name, input, output) => {
    renderSearchEntry(<WidgetHeader id="0" {...input}/>)
    await expectFilterTitle(
      output.quantity,
      output.label,
      output.description,
      output.unit,
      input.disableUnit
    )
  }
)
