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
import { readArchive, renderNoAPI } from '../../conftest.spec'
import {
  expectComposition,
  expectLatticeParameters,
  expectStructure,
  expectSymmetry,
  expectNoStructure,
  expectNoSymmetry,
  expectNoLatticeParameters
} from '../conftest.spec'
import MaterialCard from './MaterialCard'
import { cloneDeep } from 'lodash'

let archive
let index
let properties

beforeAll(async () => {
  ({archive, index, properties} = await readArchive('../../../tests/states/archives/dft.json'))
})

test('correctly renders entry with all material information', async () => {
  renderNoAPI(<MaterialCard index={index} properties={properties} archive={archive}/>)
  expectComposition(index)
  expectSymmetry(index)
  expectLatticeParameters(index)
  await expectStructure(index)
})

test('correctly renders material without symmetry information', async () => {
  const indexNew = cloneDeep(index)
  delete indexNew.results.material.symmetry
  renderNoAPI(<MaterialCard index={indexNew} properties={properties} archive={archive}/>)
  expectComposition(indexNew)
  expectNoSymmetry(indexNew)
  expectLatticeParameters(indexNew)
  await expectStructure(indexNew)
})

test('correctly renders material without lattice information', async () => {
  const indexNew = cloneDeep(index)
  delete indexNew.results.properties.structures.structure_original.lattice_parameters
  delete indexNew.results.properties.structures.structure_conventional.lattice_parameters
  delete indexNew.results.properties.structures.structure_primitive.lattice_parameters
  renderNoAPI(<MaterialCard index={indexNew} properties={properties} archive={archive}/>)
  expectComposition(indexNew)
  expectSymmetry(indexNew)
  expectNoLatticeParameters(indexNew)
  await expectStructure(indexNew)
})

test('correctly renders material without any structure information', async () => {
  const indexNew = cloneDeep(index)
  delete indexNew.results.properties.structures
  renderNoAPI(<MaterialCard index={indexNew} properties={properties} archive={archive}/>)
  expectComposition(indexNew)
  expectSymmetry(indexNew)
  expectNoLatticeParameters(indexNew)
  expectNoStructure(indexNew)
})
