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
import { renderSearchEntry, expectFilterTitle } from './conftest.spec'
import FilterTitle from './FilterTitle'

test.each([
  [
    'default',
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective'
    },
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      label: 'U effective',
      description: 'Value of the effective U parameter (u - j).',
      unit: 'eV'
    }
  ],
  [
    'custom label',
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      label: 'Custom Label'
    },
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      label: 'Custom Label',
      unit: 'eV'
    }
  ],
  [
    'custom description',
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      description: 'This is a custom description.'
    },
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      label: 'U effective',
      description: 'This is a custom description.',
      unit: 'eV'
    }
  ],
  [
    'custom all',
    {
      label: 'Custom Quantity',
      description: 'This is a completely custom quantity.',
      disableUnit: true
    },
    {
      label: 'Custom Quantity',
      description: 'This is a completely custom quantity.'
    }
  ],
  [
    'disable unit',
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      disableUnit: true
    },
    {
      quantity: 'results.method.simulation.dft.hubbard_kanamori_model.u_effective',
      label: 'U effective',
      unit: undefined
    }
  ]
])('%s', async (name, input, output) => {
    renderSearchEntry(<FilterTitle {...input}/>)
    await expectFilterTitle(
      output.quantity,
      output.label,
      output.description,
      output.unit,
      input.disableUnit
    )
  }
)
