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
import {TestAdaptor} from "./conftest.spec"

async function createMetainfo(data, parent, url = systemMetainfoUrl) {
  data._url = url
  data._metainfo = new Metainfo(parent, data, null, {}, {})
  return await data._metainfo._result
}

const mockPackage = ({
  packages: [
    {
      name: 'testPackage',
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
              m_annotations: {eln: [{component: "NumberEditQuantity"}]},
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

test.each([
  ['non-editable, no display unit', {value1: 7.5}, false, '7.5·10+10'],
  ['editable, no display unit', {value2: 7.5}, true, '75000000000']
])('Test editable/uneditable sections: %s', async (name, section, editable, expected) => {
  const metainfo = await createMetainfo(mockPackage)
  const defsByName = await metainfo.getDefsByName()
  const def = defsByName.TestSection[0]
  const adaptor = new TestAdaptor('', 'Data')

  render(
    <laneContext.Provider value={{next: {}, adaptor: adaptor}}>
      <Section
        section={section}
        def={def}
        sectionIsInEln={editable}
        sectionIsEditable={editable}
      />
    </laneContext.Provider>
  )
  if (editable) {
    const numberFieldValue = screen.getByTestId('number-edit-quantity-value')
    const numberFieldValueInput = within(numberFieldValue).getByRole('textbox')
    await waitFor(() => expect(numberFieldValueInput.value).toEqual(expected))
  } else {
    const numberWithUnit = document.querySelector('[data-testid="scientific-number-with-unit"]')
    const span = numberWithUnit.querySelector('span')
    expect(span.textContent).toContain(expected)
  }
})
