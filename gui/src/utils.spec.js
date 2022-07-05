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

import { authorList, nameList, titleCase, parseNomadUrl, normalizeNomadUrl, refType } from './utils'
import { apiBase } from './config'

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

test.each([
  ['empty value', ''],
  ['bad type, expected string, got object', {}],
  ['bad type, expected string, got number', 7],
  ['absolute nomad installation url does not contain "/api"', 'https://my.nomad.oasis.com/prod'],
  ['expected "/uploads/<uploadId>" in absolute url', 'https://my.nomad.oasis.com/prod/api/silly/continuation#more/silly'],
  ['bad start sequence', 'gibberish'],
  ['expected "/uploads/<uploadId>" in url', '../uploads'],
  ['expected "/upload" or "/uploads" in url', '../silly/ref'],
  ['archive" must be followed by entry id or "mainfile"', '../upload/archive'],
  ['"mainfile" must be followed by a mainfile path', 'https://my.nomad.oasis.com/prod/api/uploads/SomeUploadID/archive/mainfile#arch/path'],
  ['unexpected path element after entry id', '../upload/archive/SomeEntryID/silly'],
  ['expected "raw" or "archive" after upload ref', '../upload/silly'],
  ['Unexpected "#" without entry reference', '../uploads/SomeUploadId/raw#silly/arch/path']
])('parseNomadUrl fails: %s', (errorSubstring, url) => {
  expect(() => {
    parseNomadUrl(url)
  }).toThrow(errorSubstring)
})

test.each([
  ['https://my.nomad.oasis.com/prod/api', {
    originalUrlRelativeTo: null,
    type: refType.installation,
    installationUrl: 'https://my.nomad.oasis.com/prod/api',
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: undefined
  }],
  [`${apiBase}/uploads/SomeUploadID`, {
    originalUrlRelativeTo: null,
    type: refType.upload,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: undefined
  }],
  [`${apiBase}/uploads/SomeUploadID/raw`, {
    originalUrlRelativeTo: null,
    type: refType.upload,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: ''
  }],
  [`${apiBase}/uploads/SomeUploadID/raw/some/path#arch/path`, {
    originalUrlRelativeTo: null,
    type: refType.archive,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: 'arch/path'
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/SomeArchID`, {
    originalUrlRelativeTo: null,
    type: refType.archive,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: undefined
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/SomeArchID#arch/path`, {
    originalUrlRelativeTo: null,
    type: refType.archive,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: 'arch/path'
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path`, {
    originalUrlRelativeTo: null,
    type: refType.archive,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: undefined
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path#arch/path`, {
    originalUrlRelativeTo: null,
    type: refType.archive,
    installationUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: 'arch/path'
  }],
  [`../uploads/SomeUploadID`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.upload,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: undefined
  }],
  [`../uploads/SomeUploadID/raw/some/path`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.upload,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: 'some/path'
  }],
  [`../uploads/SomeUploadID/raw/some/path#definitions/some/schema/path`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.customSchema,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: 'definitions/some/schema/path'
  }],
  [`../uploads/SomeUploadID/archive/SomeArchID`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: undefined
  }],
  [`../uploads/SomeUploadID/archive/SomeArchID#arch/path`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: 'arch/path'
  }],
  [`../uploads/SomeUploadID/archive/mainfile/some/path`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: undefined
  }],
  [`../uploads/SomeUploadID/archive/mainfile/some/path#definitions/some/schema/path`, {
    originalUrlRelativeTo: refType.installation,
    type: refType.customSchema,
    installationUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: 'definitions/some/schema/path'
  }],
  [`../upload`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.upload,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: undefined
  }],
  [`../upload/raw/some/path`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.upload,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: 'some/path'
  }],
  [`../upload/raw/some/path#arch/path`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: 'some/path',
    path: 'arch/path'
  }],
  [`../upload/archive/SomeArchID`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: undefined
  }],
  [`../upload/archive/SomeArchID#definitions/path`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.customSchema,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: 'definitions/path'
  }],
  [`../upload/archive/mainfile/some/path`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: 'some/path',
    path: undefined
  }],
  [`../upload/archive/mainfile/some/path#arch/path`, {
    originalUrlRelativeTo: refType.upload,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: 'some/path',
    path: 'arch/path'
  }],
  [`#arch/path`, {
    originalUrlRelativeTo: refType.archive,
    type: refType.archive,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: 'arch/path'
  }],
  [`#definitions/def/path`, {
    originalUrlRelativeTo: refType.archive,
    type: refType.customSchema,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: 'definitions/def/path'
  }],
  [`nomad.datamodel.some.path`, {
    originalUrlRelativeTo: null,
    type: refType.globalSchema,
    installationUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: 'nomad.datamodel.some.path'
  }]
])('parseNomadUrl: %s', (url, expectedResult) => {
  const result = parseNomadUrl(url)
  expect(result.originalUrl).toBe(url)
  for (const key of ['originalUrlRelativeTo', 'type', 'installationUrl', 'uploadId', 'entryId', 'mainfile', 'path']) {
    expect(key in result).toBe(true)
    expect(result[key]).toBe(expectedResult[key])
  }
})

test.each([
  ['a baseUrl is required', '../upload', undefined],
  ['baseUrl is not absolute', '../upload', '../uploads/SomeUploadID'],
  ['missing information about uploadId', '../upload/raw/some/path', apiBase],
  ['missing information about entryId', '#some/path', `${apiBase}/uploads/SomeUploadID`]
])('normalizeNomadUrl fails: %s', (errorSubstring, url, baseUrl) => {
  expect(() => {
    normalizeNomadUrl(url, baseUrl)
  }).toThrow(errorSubstring)
})

test.each([
  ['../uploads/SomeUploadID/archive/SomeArchID#x/y/z', apiBase, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeArchID#x/y/z`
  }],
  ['../uploads/SomeUploadID/raw/some/path', `${apiBase}/uploads/SomeOtherUploadID/raw/other/path`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/raw/some/path`
  }],
  ['../uploads/SomeUploadID/raw/some/path#arch/path', `${apiBase}/uploads/SomeOtherUploadID/raw/other/path`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/TbJz7EfLcUdPBJ_iSAXrm5cy7G1v#arch/path`
  }],
  ['../uploads/SomeUploadID/archive/SomeArchID#x/y/z', `${apiBase}/uploads/SomeOtherUploadID/archive/SomeOtherArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeArchID#x/y/z`
  }],
  ['../upload/raw/some/path', `${apiBase}/uploads/SomeUploadID/archive/SomeArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/raw/some/path`
  }],
  ['../upload/raw/some/path#arch/path', `${apiBase}/uploads/SomeUploadID/archive/SomeArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/TbJz7EfLcUdPBJ_iSAXrm5cy7G1v#arch/path`
  }],
  ['../upload/raw/some/path', `https://other.nomad.com/nomd/api/uploads/SomeUploadID/archive/SomeArchID`, {
    normalizedUrl: `https://other.nomad.com/nomd/api/uploads/SomeUploadID/raw/some/path`,
    isExternal: true
  }],
  ['../upload/archive/SomeArchID#definitions/x/y/z', `${apiBase}/uploads/SomeUploadID/archive/SomeOtherArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeArchID#definitions/x/y/z`
  }],
  ['#x/y/z', `${apiBase}/uploads/SomeUploadID/archive/SomeOtherArchID#a/b/c`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeOtherArchID#x/y/z`
  }],
  ['#x/y/z', `${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path#a/b/c`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/TbJz7EfLcUdPBJ_iSAXrm5cy7G1v#x/y/z`
  }],
  ['nomad.datamodel.some.path', `${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path#a/b/c`, {
    normalizedUrl: `nomad.datamodel.some.path`
  }],
  ['https://other.nomad.com/nomd/api/uploads/SomeUploadID/raw/x/y/z', `${apiBase}/uploads/SomeOtherUploadID/archive/mainfile/some/path#a/b/c`, {
    normalizedUrl: `https://other.nomad.com/nomd/api/uploads/SomeUploadID/raw/x/y/z`,
    isExternal: true
  }]
])('normalizeNomadUrl url = "%s", baseUrl = "%s"', (url, baseUrl, expectedResult) => {
  const result = normalizeNomadUrl(url, baseUrl)
  expect(result.normalizedUrl).toBe(expectedResult.normalizedUrl)
  expect(result.isExternal).toBe(expectedResult.isExternal || false)
})
