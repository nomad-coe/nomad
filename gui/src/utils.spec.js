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

import { authorList, nameList, titleCase } from './utils'

describe('titleCase', () => {
  it('runs on empty strings', () => {
    expect(titleCase('')).toBe('')
  })
  it('upcases each word', () => {
    expect(titleCase('a')).toBe('A')
    expect(titleCase('ab')).toBe('Ab')
    expect(titleCase('a b')).toBe('A B')
    expect(titleCase('aa bb')).toBe('Aa Bb')
    expect(titleCase('aa_bb')).toBe('Aa Bb')
  })
})

const longUserList = [{name: 'a'}, {name: 'b'}, {name: 'c'}, {name: 'd'}]
describe('nameList', () => {
  it('abbreviates after 2', () => {
    expect(nameList(longUserList)).toBe('A, D et al')
  })
  it('does not abbreviate with expanded', () => {
    expect(nameList(longUserList, true)).toBe('A, B, C, D')
  })
})

describe('authorList', () => {
  it('runs on empty entry', () => {
    expect(authorList(null)).toBe('')
  })
  it('returns only external db', () => {
    expect(authorList({external_db: 'AFLOW', authors: longUserList})).toBe('AFLOW')
  })
  it('returns author list', () => {
    expect(authorList({authors: longUserList})).toBe('A, D et al')
  })
  it('expands external db with short author list', () => {
    expect(authorList({external_db: 'AFLOW', authors: longUserList}, true)).toBe('AFLOW (A, D et al)')
  })
  it('expands author list', () => {
    expect(authorList({authors: longUserList}, true)).toBe('A, B, C, D')
  })
})
