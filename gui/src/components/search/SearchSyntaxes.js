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

// Contains definitions for the different search bar syntaxes. This mapping is
// used for performing the parsing, creating an info dialog and for
// including/excluding certain syntaxes with the app config.

const opMap = {
  '<=': 'lte',
  '>=': 'gte',
  '>': 'gt',
  '<': 'lt'
}
const opMapReverse = {
  '<=': 'gte',
  '>=': 'lte',
  '>': 'lt',
  '<': 'gt'
}

const reString = '[^\\s=<>](?:[^=<>]*[^\\s=<>])?'
const reOp = '(?:<|>)=?'
const reEquality = new RegExp(`^(${reString})\\s*=\\s*(${reString})$`)
const reLte = RegExp(`^(${reString})\\s*(${reOp})\\s*(${reString})$`)
const reLteGte = new RegExp(`^(${reString})\\s*(${reOp})\\s*(${reString})\\s*(${reOp})\\s*(${reString})$`)
const reExists = new RegExp(`^(${reString})\\s*=\\s*\\*$`)

export const SearchSyntaxes = {
  existence: {
    regex: reExists,
    parse: (input) => {
      const exists = input.match(reExists)
      return {
        inputNormalized: `${exists[1]} = *`,
        target: 'quantities',
        value: exists[1]
      }
    },
    readme: 'Used to query for the existence of a specific metainfo field in the data.',
    examples: ['authors.name = *'],
    label: 'Value existence',
    labelShort: 'existence'
  },
  equality: {
    regex: reEquality,
    parse: (input) => {
      const equals = input.match(reEquality)
      const target = equals[1]
      return {
        inputNormalized: `${equals[1]} = ${equals[2]}`,
        target,
        value: equals[2]
      }
    },
    readme: 'Used to query for a specific value with exact match.',
    examples: ['author.name = John Doe'],
    label: 'Value equality',
    labelShort: 'equality'
  },
  range_bounded: {
    regex: reLteGte,
    parse: (input) => {
      const ltegteSandwich = input.match(reLteGte)
      const a = ltegteSandwich[1]
      const op1 = ltegteSandwich[2]
      const b = ltegteSandwich[3]
      const op2 = ltegteSandwich[4]
      const c = ltegteSandwich[5]
      const target = b
      const value = {
        [opMapReverse[op1]]: a,
        [opMap[op2]]: c
      }
      return {
        inputNormalized: `${ltegteSandwich[1]} ${ltegteSandwich[2]} ${ltegteSandwich[3]} ${ltegteSandwich[4]} ${ltegteSandwich[5]}`,
        target,
        value
      }
    },
    readme: 'Queries values that are between two numerical limits, inclusive or exclusive. You can also specify a unit if the target has a certain dimensionality.',
    examples: [
      '0 <= results.material.n_elements <= 2',
      '0 eV < results.properties.geometry_optimization.final_energy_difference < 1e-3 eV'
    ],
    label: 'Range query (bounded)',
    labelShort: 'range'
  },
  range_half_bounded: {
    regex: reLte,
    parse: (input) => {
      const ltegte = input.match(reLte)
      const a = ltegte[1]
      const op = ltegte[2]
      const b = ltegte[3]
      const aNumber = /\d/.test(a)
      const target = aNumber ? b : a
      const query = aNumber ? a : b
      const value = {[opMap[op]]: query}
return {
        inputNormalized: `${ltegte[1]} ${ltegte[2]} ${ltegte[3]}`,
        target,
        value
      }
    },
    readme: 'Queries values that are above/below a numerical limit, inclusive or exclusive. You can also specify a unit if the target has a certain dimensionality.',
    examples: [
      'results.material.n_elements > 2',
      'results.properties.geometry_optimization.final_energy_difference <= 1e-3 eV'
    ],
    label: 'Range query (half-bounded)',
    labelShort: 'range'
  },
  free_text: {
    regex: /[^]*/,
    parse: (input) => {
      return {
        inputNormalized: input,
        target: 'text_search_contents',
        value: input
      }
    },
    readme: 'Free-form text query across keywords associated with an entry. Requires that a set of text search contents has been filled in the entry.',
    examples: ['type anything there'],
    label: 'Free-text query',
    labelShort: 'free-text'
  }
}
