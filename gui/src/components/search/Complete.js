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

import searchQuantities from '../../searchQuantities'

// Mapping from full name to abbreviation, if abbreviation is available
export const quantityAbbreviations = new Map()

// Mapping from abbreviation to full name, if abbreviation is available
export const quantityFullnames = new Map()

const nonUniques = new Set()
const quantityNames = Object.keys(searchQuantities).filter(q => q.startsWith('results'))
const abbreviationList = quantityNames.map(name => name.split('.').pop())
const abbreviationSet = new Set()
for (let abbr of abbreviationList) {
  if (abbreviationSet.has(abbr)) {
    nonUniques.add(abbr)
  }
  abbreviationSet.add(abbr)
}
for (let i = 0; i < quantityNames.length; ++i) {
  const abbr = abbreviationList[i]
  const name = quantityNames[i]
  const unique = !nonUniques.has(abbr)
  if (unique) {
    quantityAbbreviations.set(name, abbr)
    quantityFullnames.set(abbr, name)
  }
}

export const opMap = {
  '<=': 'lte',
  '>=': 'gte',
  '>': 'gt',
  '<': 'lt'
}
export const opMapReverse = {
  '<=': 'gte',
  '>=': 'lte',
  '>': 'lt',
  '<': 'gt'
}
