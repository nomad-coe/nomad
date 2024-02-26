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
import {render, screen} from "../conftest.spec"
import {within} from "@testing-library/dom"
import {waitFor} from "@testing-library/react"
import {Section} from "./ArchiveBrowser"
import {Metainfo} from "./metainfo"
import {systemMetainfoUrl} from "../../utils"
import {laneContext} from './Browser'
import {TestAdaptor} from "./Browser.spec"

async function createMetainfo(data, parent, url = systemMetainfoUrl) {
  data._url = url
  data._metainfo = new Metainfo(parent, data, null, {}, {})
  return await data._metainfo._result
}

const mockPackage = ({
  packages: [
    {
      name: 'testPackage',
      m_annotations: {
        display: [{
          unit_system: 'AU'
        }]
      },
      section_definitions: [
        {
          name: 'TestSection',
          quantities: [
            {
              name: 'value1',
              type: { type_kind: 'python', type_data: 'float' },
              unit: 'meter',
              m_parent_sub_section: 'quantities'
            },
            {
              name: 'value2',
              m_annotations: {
                eln: [
                  {
                    component: "NumberEditQuantity"
                  }
                ]
              },
              type: { type_kind: 'python', type_data: 'float' },
              unit: 'meter',
              m_parent_sub_section: 'quantities'
            }
          ]
        },
        {
          name: 'TestSection',
          m_annotations: {
            display: [{
              unit_system: 'SI'
            }]
          },
          quantities: [
            {
              name: 'value3',
              type: { type_kind: 'python', type_data: 'float' },
              unit: 'meter',
              m_parent_sub_section: 'quantities'
            },
            {
              name: 'value4',
              m_annotations: {
                eln: [
                  {
                    component: "NumberEditQuantity"
                  }
                ]
              },
              type: { type_kind: 'python', type_data: 'float' },
              unit: 'meter',
              m_parent_sub_section: 'quantities'
            }
          ]
        }
      ]
    }
  ]
})

describe('Test interaction between unit menu and quantities', () => {
  it('Package with display unit_system', async () => {
    const metainfo = await createMetainfo(mockPackage)
    const defsByName = await metainfo.getDefsByName()
    const def = defsByName.TestSection[0]
    const adaptor = new TestAdaptor('', 'Data')

    render(
      <laneContext.Provider value={{next: {}, adaptor: adaptor}}>
        <Section
          section={{value1: 7.5}}
          def={def}
          sectionIsInEln={false}
          sectionIsEditable={false}
        />
      </laneContext.Provider>
    )
    // should be rendered in package default unit system
    const numberWithUnit = document.querySelector('[data-testid="scientific-number-with-unit"]')
    const span = numberWithUnit.querySelector('span')
    expect(span.textContent).toContain('1.41729459')
    expect(span.textContent).toContain('·10')
    const sup = document.querySelector('sup')
    expect(sup.textContent).toContain('+11')
  })

  it('Editable package with display unit_system', async () => {
    const metainfo = await createMetainfo(mockPackage)
    const defsByName = await metainfo.getDefsByName()
    const def = defsByName.TestSection[0]
    const adaptor = new TestAdaptor('', 'Data')

    render(
      <laneContext.Provider value={{next: {}, adaptor: adaptor}}>
        <Section
          section={{value2: 8.5}}
          def={def}
          sectionIsInEln={true}
          sectionIsEditable={true}
        />
      </laneContext.Provider>
    )
    // should be rendered in package default unit system
    const numberFieldValue = screen.getByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    // should be rendered in default unit system
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('160626720592.894'))
  })

  it('Section with display unit_system', async () => {
    const metainfo = await createMetainfo(mockPackage)
    const defsByName = await metainfo.getDefsByName()
    const def = defsByName.TestSection[1]
    const adaptor = new TestAdaptor('', 'Data')

    render(
      <laneContext.Provider value={{next: {}, adaptor: adaptor}}>
        <Section
          section={{value3: 7.5}}
          def={def}
          sectionIsInEln={false}
          sectionIsEditable={false}
        />
      </laneContext.Provider>
    )
    // Should be rendered in section default unit system, The determined package default unit system should be overridden.
    const numberWithUnit = document.querySelector('[data-testid="scientific-number-with-unit"]')
    expect(numberWithUnit.textContent).toContain('7.50000')
  })

  it('Editable Section with display unit_system', async () => {
    const metainfo = await createMetainfo(mockPackage)
    const defsByName = await metainfo.getDefsByName()
    const def = defsByName.TestSection[1]
    const adaptor = new TestAdaptor('', 'Data')

    render(
      <laneContext.Provider value={{next: {}, adaptor: adaptor}}>
        <Section
          section={{value4: 8.5}}
          def={def}
          sectionIsInEln={true}
          sectionIsEditable={true}
        />
      </laneContext.Provider>
    )
    // Should be rendered in section default unit system, The determined package default unit system should be overridden.
    const numberFieldValue = screen.getByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    await waitFor(() => expect(numberFieldValueInput.value).toEqual('8.5'))
  })
})
