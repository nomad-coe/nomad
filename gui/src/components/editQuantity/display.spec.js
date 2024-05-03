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
import {editQuantityComponents} from "./EditQuantity"
import {PackageMDef, QuantityMDef} from "../archive/metainfo"
import {StringEditQuantity} from "./StringEditQuantity"
import ListEditQuantity from "./ListEditQuantity"

test.each([
  [
    'NumberEditQuantity without custom label',
    'NumberEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'NumberEditQuantity with custom label',
    'NumberEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'StringEditQuantity without custom label',
    'StringEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'StringEditQuantity with custom label',
    'StringEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'URLEditQuantity without custom label',
    'URLEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'URLEditQuantity with custom label',
    'URLEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'EnumEditQuantity without custom label',
    'EnumEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'EnumEditQuantity with custom label',
    'EnumEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'SelectEnumEditQuantity without custom label',
    'SelectEnumEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'SelectEnumEditQuantity with custom label',
    'SelectEnumEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'RadioEnumEditQuantity without custom label',
    'RadioEnumEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'RadioEnumEditQuantity with custom label',
    'RadioEnumEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'AutocompleteEditQuantity without custom label',
    'AutocompleteEditQuantity',
    {type: {type_kind: 'Enum', type_data: []}},
    undefined,
    'my value'
  ],
  [
    'AutocompleteEditQuantity with custom label',
    'AutocompleteEditQuantity',
    {type: {type_kind: 'Enum', type_data: []}},
    'My custom value',
    'My custom value'
  ],
  [
    'BoolEditQuantity without custom label',
    'BoolEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'BoolEditQuantity with custom label',
    'BoolEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'FileEditQuantity without custom label',
    'FileEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'FileEditQuantity with custom label',
    'FileEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'DateTimeEditQuantity without custom label',
    'DateTimeEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'DateTimeEditQuantity with custom label',
    'DateTimeEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'RichTextEditQuantity without custom label',
    'RichTextEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'RichTextEditQuantity with custom label',
    'RichTextEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ],
  [
    'ReferenceEditQuantity without custom label',
    'ReferenceEditQuantity',
    {type: {_referencedDefinition: {m_def: QuantityMDef, _section: {m_def: PackageMDef, _allBaseSections: [], _qualifiedName: 'a', _pkgParentData: {_metainfo: {_parsedUrl: 'a'}}}}}},
    undefined,
    'my value'
  ],
  [
    'ReferenceEditQuantity with custom label',
    'ReferenceEditQuantity',
    {type: {_referencedDefinition: {m_def: QuantityMDef, _section: {m_def: PackageMDef, _allBaseSections: [], _qualifiedName: 'a', _pkgParentData: {_metainfo: {_parsedUrl: 'a'}}}}}},
    'My custom value',
    'My custom value'
  ],
  [
    'AuthorEditQuantity without custom label',
    'AuthorEditQuantity',
    undefined,
    undefined,
    'my value'
  ],
  [
    'AuthorEditQuantity with custom label',
    'AuthorEditQuantity',
    undefined,
    'My custom value',
    'My custom value'
  ]
])('Test label of %s', async (name, component, def, customLabel, expected) => {
  render(React.createElement(
    editQuantityComponents[component], {
      quantityDef: {
        name: 'my_value',
        m_annotations: {display: [{label: customLabel}]},
        ...def
      },
      type: {type_kind: 'python', type_data: 'string'},
      value: undefined
    }
    ))
  screen.getByText(expected)
})

test.each([
  [
    'without custom label',
    undefined,
    'my value'
  ],
  [
    'with custom label',
    'My custom value',
    'My custom value'
  ]
])('Test label of ListEditQuantity %s', async (name, customLabel, expected) => {
  render(
    <ListEditQuantity
      component={StringEditQuantity}
      quantityDef={{
        name: 'my_value',
        m_annotations: {display: [{label: customLabel}]}
      }}
      type={{type_kind: 'python', type_data: 'string'}}
    />)
  const labels = screen.getAllByText(expected)
  expect(labels.length).toBe(2)
})
