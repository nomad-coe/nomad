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

import { createMetainfo } from './metainfo'

test('metainfo initializes correctly', async () => {
  const data = ({
    packages: [
      {
        name: 'testPackage',
        section_definitions: [
          {
            name: 'TestSection',
            base_sections: [
              '/packages/0/section_definitions/1'
            ],
            quantities: [
              {
                name: 'testQuantity',
                type: { type_kind: 'python', type_data: 'str' }
              }
            ]
          },
          {
            name: 'TestBaseSection',
            quantities: [
              {
                name: 'baseSectionQuantity',
                type: { type_kind: 'python', type_data: 'str' }
              }
            ]
          }
        ]
      }
    ]
  })

  const metainfo = await createMetainfo(data)

  expect(typeof metainfo).toBe('object')
  const defsByName = await metainfo.getDefsByName()
  expect(defsByName.TestSection.length).toBe(1)
  expect(defsByName.TestSection[0]._allBaseSections.length).toBe(1)
  expect(defsByName.TestSection[0]._allProperties.length).toBe(2)
  expect(defsByName.TestSection[0]._properties.testQuantity.name).toBe('testQuantity')
  expect(defsByName.TestSection[0]._properties.baseSectionQuantity.name).toBe('baseSectionQuantity')
})

test('parent metainfo definitions are resolved', async () => {
  const parentData = ({
    packages: [
      {
        name: 'testPackage',
        section_definitions: [
          {
            name: 'TestSection',
            quantities: [
              {
                name: 'testQuantity',
                type: { type_kind: 'python', type_data: 'str' }
              }
            ]
          }
        ]
      }
    ]
  })

  const childData = ({
    definitions: {
      m_def: 'nomad.metainfo.metainfo.Package',
      section_definitions: [
        {
          name: 'DerivedTestSection',
          base_sections: [
            'testPackage.TestSection'
          ]
        }
      ]
    }
  })

  const parentMetainfo = await createMetainfo(parentData)
  const childMetainfo = await createMetainfo(childData, parentMetainfo)

  const defsByName = await childMetainfo.getDefsByName()
  expect(defsByName.DerivedTestSection[0]._properties.testQuantity.name).toBe('testQuantity')
})
