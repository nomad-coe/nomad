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
import { range } from 'lodash'
import { TestAdaptor } from './conftest.spec'
import { screen, renderNoAPI, render } from '../conftest.spec'
import { expectPagination } from '../visualization/conftest.spec'
import { PropertyValuesList, Section } from './ArchiveBrowser'
import { laneContext } from './Browser'
import {waitFor} from "@testing-library/dom"
import {Metainfo} from './metainfo'
import {systemMetainfoUrl} from '../../utils'

test.each([
  [15, 10, 5],
  [12, 10, 5]
])('test subsection with no pagination, items: %s, top: %s, bottom:%s', async (nItems, nTop, nBottom) => {
  const indices = range(nItems)
  const label = "subsection"
  const values = indices.map(i => (null))

  renderNoAPI(
    <laneContext.Provider value={{next: {}}}>
      <PropertyValuesList itemKey={label} values={values} nTop={nTop} nBottom={nBottom} open={true} />
    </laneContext.Provider>
  )

  // Expect to find all items
  for (const i of indices) {
    screen.getByText(`${i}`)
  }
  // Pagination component should not be visible
  await expectPagination(false, false, false)
})

test.each([
  [30, 10, 5],
  [16, 10, 5],
  [30, 10, 0]
])('test subsection with pagination, items: %s, top: %s, bottom:%s', async (nItems, nTop, nBottom) => {
  const indices = range(nItems)
  const label = "subsection"
  const values = indices.map(i => (null))

  renderNoAPI(
    <laneContext.Provider value={{next: {}}}>
      <PropertyValuesList itemKey={label} values={values} nTop={nTop} nBottom={nBottom} open={true} />
    </laneContext.Provider>
  )

  // Expect to find top and bottom items
  for (const i of range(nTop)) {
    screen.getByText(`${i}`)
  }
  for (const i of range(nBottom)) {
    screen.getByText(`${nItems - i - 1}`)
  }
  // Both pagination components should be visible
  const downPagination = screen.getByTestId('propertyvalueslist-pagination-down')
  await expectPagination(true, false, false, downPagination)
  if (nBottom > 0) {
    const upPagination = screen.getByTestId('propertyvalueslist-pagination-up')
    await expectPagination(true, false, false, upPagination)
  } else {
    await waitFor(() => expect(screen.queryByTestId('propertyvalueslist-pagination-up')).not.toBeInTheDocument())
  }
})

async function createMetainfo(data, parent, url = systemMetainfoUrl) {
  data._url = url
  data._metainfo = new Metainfo(parent, data, null, {}, {})
  return await data._metainfo._result
}

const mockPackage = (properties) => {
  const stringEditQuantity = (name) => (
    {
      name: name,
      default: name,
      m_annotations: {
        eln: [
          {
            component: "StringEditQuantity"
          }
        ]
      },
      m_parent_sub_section: 'quantities'
    }
  )

  return {
  packages: [
    {
      name: 'testPackage',
      section_definitions: [
        {
          name: 'TestSection',
          m_annotations: {
            eln: [
              {
                properties: properties
              }
            ]
          },
          quantities: ['value1', 'value2', 'value3'].map(stringEditQuantity)
        }
      ]
    }
  ]
}
}

test.each([
  [
    'visible included',
    {
      visible: {
        include: ['value2']
      }
    },
    ['value2']
  ],
  [
    'visible excluded',
    {
      visible: {
        exclude: ['value2']
      }
    },
    ['value1', 'value3']
  ],
  [
    'visible included and excluded',
    {
      visible: {
        include: ['value1', 'value3'],
        exclude: ['value1']
      }
    },
    ['value3']
  ]
])('test eln properties: %s', async (name, properties, expected) => {
  const metainfo = await createMetainfo(mockPackage(properties))
  const defsByName = await metainfo.getDefsByName()
  const def = defsByName.TestSection[0]
  const adaptor = new TestAdaptor('', 'Data')

  render(
    <laneContext.Provider value={{next: {}, adaptor: adaptor}}>
      <Section
        section={{}}
        def={def}
        sectionIsInEln={true}
        sectionIsEditable={true}
      />
    </laneContext.Provider>
  )

  const numberFields = screen.getAllByRole('textbox')
  expect(numberFields).toHaveLength(expected.length)
  for (const [index, value] of expected.entries()) {
    await waitFor(() => expect(numberFields[index].value).toEqual(value))
  }
})
