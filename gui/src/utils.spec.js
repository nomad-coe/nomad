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

import {
  authorList,
  nameList,
  titleCase,
  parseNomadUrl,
  resolveNomadUrl,
  normalizeNomadUrl,
  refType,
  refRelativeTo,
  parseJMESPath
} from './utils'
import { apiBase, urlAbs } from './config'
import { isEqual } from 'lodash'

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
  ['absolute nomad deployment url does not contain "/api"', 'https://my.nomad.oasis.com/prod'],
  ['expected "/uploads/<uploadId>" in absolute url', 'https://my.nomad.oasis.com/prod/api/silly/continuation#/more/silly'],
  ['bad url', 'a.b.0c'],
  ['expected "/uploads/<uploadId>" in url', '../uploads'],
  ['expected "/upload" or "/uploads" in url', '../silly/ref'],
  ['archive" must be followed by entry id or "mainfile"', '../upload/archive'],
  ['"mainfile" must be followed by a mainfile path', 'https://my.nomad.oasis.com/prod/api/uploads/SomeUploadID/archive/mainfile#/arch/path'],
  ['unexpected path element after entry id', '../upload/archive/SomeEntryID/silly'],
  ['expected "raw" or "archive" after upload ref', '../upload/silly'],
  ['Unexpected dataPath without entry reference', '../uploads/SomeUploadId/raw#/silly/arch/path'],
  ['versionHash can only be specified for metainfo urls.', '../uploads/SomeUploadID/archive/SomeArchID#/arch/path@SomeHash'],
  ['bad versionHash provided', '../uploads/SomeUploadID/archive/SomeArchID#/definitions/some/path@*']

])('parseNomadUrl fails: %s', (errorSubstring, url) => {
  expect(() => {
    parseNomadUrl(url)
  }).toThrow(errorSubstring)
})

test.each([
  ['https://my.nomad.oasis.com/prod/api', {
    relativeTo: null,
    type: refType.deployment,
    deploymentUrl: 'https://my.nomad.oasis.com/prod/api',
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: true
  }],
  [`${apiBase}/uploads/SomeUploadID`, {
    relativeTo: null,
    type: refType.upload,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`${apiBase}/uploads/SomeUploadID/raw`, {
    relativeTo: null,
    type: refType.upload,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: '',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`${apiBase}/uploads/SomeUploadID/raw/some/path#/arch/path`, {
    relativeTo: null,
    type: refType.archive,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/SomeArchID`, {
    relativeTo: null,
    type: refType.archive,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/SomeArchID#arch/path`, {
    relativeTo: null,
    type: refType.archive,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path`, {
    relativeTo: null,
    type: refType.archive,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path#/arch//path`, {
    relativeTo: null,
    type: refType.archive,
    deploymentUrl: apiBase,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: true,
    isExternal: false
  }],
  [`../uploads/SomeUploadID`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.upload,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/raw/some/path`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.upload,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: undefined,
    mainfile: undefined,
    path: 'some/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/raw/some/path#/definitions/some/schema/path`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: '/definitions/some/schema/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/raw/some/path#definitions/some/schema/path@SomeVersionHash`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: '/definitions/some/schema/path',
    qualifiedName: undefined,
    versionHash: 'SomeVersionHash',
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/archive/SomeArchID`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/archive/SomeArchID#/arch/path`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/archive/mainfile/some/path`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/archive/mainfile/some/path#/definitions/some/schema/path`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: '/definitions/some/schema/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../uploads/SomeUploadID/archive/mainfile/some/path#definitions/some//schema/path@SomeVersionHash`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: 'SomeUploadID',
    entryId: 'TbJz7EfLcUdPBJ_iSAXrm5cy7G1v',
    mainfile: 'some/path',
    path: '/definitions/some/schema/path',
    qualifiedName: undefined,
    versionHash: 'SomeVersionHash',
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload`, {
    relativeTo: refRelativeTo.upload,
    type: refType.upload,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/raw/some/path`, {
    relativeTo: refRelativeTo.upload,
    type: refType.upload,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: 'some/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/raw/some/path#/arch/path`, {
    relativeTo: refRelativeTo.upload,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: 'some/path',
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/archive/SomeArchID`, {
    relativeTo: refRelativeTo.upload,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/archive/SomeArchID#/definitions/path`, {
    relativeTo: refRelativeTo.upload,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: '/definitions/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/archive/SomeArchID#/definitions/path@SomeVersionHash`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: 'SomeArchID',
    mainfile: undefined,
    path: '/definitions/path',
    qualifiedName: undefined,
    versionHash: 'SomeVersionHash',
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/archive/mainfile/some/path`, {
    relativeTo: refRelativeTo.upload,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: 'some/path',
    path: undefined,
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`../upload/archive/mainfile/some/path#/arch/path`, {
    relativeTo: refRelativeTo.upload,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: 'some/path',
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`#/arch/path`, {
    relativeTo: refRelativeTo.data,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`#/definitions/def/path`, {
    relativeTo: refRelativeTo.data,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: '/definitions/def/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`/arch/path`, {
    relativeTo: refRelativeTo.data,
    type: refType.archive,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: '/arch/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`/definitions/def/path`, {
    relativeTo: refRelativeTo.data,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: '/definitions/def/path',
    qualifiedName: undefined,
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }],
  [`nomad.datamodel.some.path`, {
    relativeTo: refRelativeTo.deployment,
    type: refType.metainfo,
    deploymentUrl: undefined,
    uploadId: undefined,
    entryId: undefined,
    mainfile: undefined,
    path: undefined,
    qualifiedName: 'nomad.datamodel.some.path',
    versionHash: undefined,
    isResolved: false,
    isExternal: undefined
  }]
])('parseNomadUrl: %s', (url, expectedResult) => {
  const result = parseNomadUrl(url)
  expect(result.url).toBe(url)
  for (const key of ['relativeTo', 'type', 'deploymentUrl', 'uploadId', 'entryId', 'mainfile', 'path']) {
    expect(key in result).toBe(true)
    expect(result[key]).toBe(expectedResult[key])
  }
})

test.each([
  ['a baseUrl is required', '../upload', undefined],
  ['unresolved baseUrl', '../upload', '../uploads/SomeUploadID'],
  ['missing information about uploadId', '../upload/raw/some/path', apiBase],
  ['missing information about entryId', '#/some/path', `${apiBase}/uploads/SomeUploadID`]
])('normalizeNomadUrl fails: %s', (errorSubstring, url, baseUrl) => {
  expect(() => {
    resolveNomadUrl(url, baseUrl)
  }).toThrow(errorSubstring)
})

test.each([
  ['../uploads/SomeUploadID/archive/SomeArchID#/x/y/z', apiBase, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeArchID#/x/y/z`
  }],
  ['../uploads/SomeUploadID/raw/some/path', `${apiBase}/uploads/SomeOtherUploadID/raw/other/path`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/raw/some/path`
  }],
  ['../uploads/SomeUploadID/raw/some/path#/arch/path', `${apiBase}/uploads/SomeOtherUploadID/raw/other/path`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/TbJz7EfLcUdPBJ_iSAXrm5cy7G1v#/arch/path`
  }],
  ['../uploads/SomeUploadID/archive/SomeArchID#/x/y/z', `${apiBase}/uploads/SomeOtherUploadID/archive/SomeOtherArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeArchID#/x/y/z`
  }],
  ['../upload/raw/some/path', `${apiBase}/uploads/SomeUploadID/archive/SomeArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/raw/some/path`
  }],
  ['../upload/raw/some/path#/arch/path', `${apiBase}/uploads/SomeUploadID/archive/SomeArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/TbJz7EfLcUdPBJ_iSAXrm5cy7G1v#/arch/path`
  }],
  ['../upload/raw/some/path', `https://other.nomad.com/nomd/api/uploads/SomeUploadID/archive/SomeArchID`, {
    normalizedUrl: `https://other.nomad.com/nomd/api/uploads/SomeUploadID/raw/some/path`,
    isExternal: true
  }],
  ['../upload/archive/SomeArchID#/definitions/x/y/z', `${apiBase}/uploads/SomeUploadID/archive/SomeOtherArchID`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeArchID#/definitions/x/y/z`
  }],
  ['#/x/y/z', `${apiBase}/uploads/SomeUploadID/archive/SomeOtherArchID#/a/b/c`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/SomeOtherArchID#/x/y/z`
  }],
  ['#/x/y/z', `${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path#/a/b/c`, {
    normalizedUrl: `${apiBase}/uploads/SomeUploadID/archive/TbJz7EfLcUdPBJ_iSAXrm5cy7G1v#/x/y/z`
  }],
  ['nomad.datamodel.some.path', `${apiBase}/uploads/SomeUploadID/archive/mainfile/some/path#/a/b/c`, {
    normalizedUrl: `nomad.datamodel.some.path`
  }],
  ['https://other.nomad.com/nomd/api/uploads/SomeUploadID/raw/x/y/z', `${apiBase}/uploads/SomeOtherUploadID/archive/mainfile/some/path#/a/b/c`, {
    normalizedUrl: `https://other.nomad.com/nomd/api/uploads/SomeUploadID/raw/x/y/z`,
    isExternal: true
  }]
])('normalizeNomadUrl url = "%s", baseUrl = "%s"', (url, baseUrl, expectedResult) => {
  const resolvedUrl = resolveNomadUrl(url, baseUrl)
  expect(resolvedUrl.isExternal).toBe(expectedResult.isExternal || false)
  expect(normalizeNomadUrl(resolvedUrl)).toBe(expectedResult.normalizedUrl)
})

test.each([
  ['already absolute', 'https://nomad-lab.eu/test', 'https://nomad-lab.eu/test', 'https://nomad-lab.eu'],
  ['root path', '/test', 'https://nomad-lab.eu/test', 'https://nomad-lab.eu'],
  ['relative path', '../test', 'https://nomad-lab.eu/test', 'https://nomad-lab.eu/folder'],
  ['protocol change', '/test', 'http://nomad-lab.eu/test', 'https://nomad-lab.eu', 'http:'],
  ['non-http protocol', 'ssh://nomad-lab.eu/test', 'ssh://nomad-lab.eu/test', 'ssh://nomad-lab.eu']
])('absolute url creation: %s ', (id, input, output, base, protocol = undefined) => {
  expect(urlAbs(input, base, protocol)).toBe(output)
})

test.each([
  [
    'simple subexpression',
    'results.material.n_elements',
    {
      quantity: 'results.material.n_elements',
      path: 'results.material.n_elements',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'index expression',
    'results.material.elements[0]',
    {
      quantity: 'results.material.elements',
      path: 'results.material.elements[0]',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'slicing ',
    'results.material.elements[0:5]',
    {
      quantity: 'results.material.elements',
      path: 'results.material.elements[0:5]',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'list projection',
    'results.properties.electronic.band_gap[*].value',
    {
      quantity: 'results.properties.electronic.band_gap.value',
      path: 'results.properties.electronic.band_gap[*].value',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'flatten projection',
    'results.properties[].electronic.band_gap[].value',
    {
      quantity: 'results.properties.electronic.band_gap.value',
      path: 'results.properties[].electronic.band_gap[].value',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'function with one argument',
    'min(results.properties.electronic.band_gap[*].value)',
    {
      quantity: 'results.properties.electronic.band_gap.value',
      path: 'min(results.properties.electronic.band_gap[*].value)',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'function with two arguments',
    'min_by(results.properties.electronic.band_gap[*], &value).type',
    {
      quantity: 'results.properties.electronic.band_gap.type',
      path: 'min_by(results.properties.electronic.band_gap[*], &value).type',
      extras: ['results.properties.electronic.band_gap.value'],
      error: undefined,
      schema: ''
    }
  ],
  [
    'filter projection',
    "results.material.topology[?label=='original'].cell.a",
    {
      quantity: 'results.material.topology.cell.a',
      path: "results.material.topology[?label=='original'].cell.a",
      extras: ['results.material.topology.label'],
      error: undefined,
      schema: ''
    }
  ],
  [
    'pipe',
    'results.properties.electronic.band_gap[*].value | min(@)',
    {
      quantity: 'results.properties.electronic.band_gap.value',
      path: 'results.properties.electronic.band_gap[*].value | min(@)',
      extras: [],
      error: undefined,
      schema: ''
    }
  ],
  [
    'schema name and dtype are handled correctly 1',
    'min_by(results.properties.electronic.band_gap[*], &value).type#MySchema#int',
    {
      quantity: 'results.properties.electronic.band_gap.type#MySchema#int',
      path: 'min_by(results.properties.electronic.band_gap[*], &value).type',
      extras: ['results.properties.electronic.band_gap.value#MySchema#int'],
      error: undefined,
      schema: '#MySchema#int'
    }
  ],
  [
    'schema name and dtype are handled correctly 2',
    'results.properties.electronic.band_gap[*].value#MySchema | min(@)',
    {
      quantity: 'results.properties.electronic.band_gap.value#MySchema',
      path: 'results.properties.electronic.band_gap[*].value | min(@)',
      extras: [],
      error: undefined,
      schema: '#MySchema'
    }
  ],
  [
    'syntax error',
    'results.material.n_elements[*',
    {
      quantity: undefined,
      path: undefined,
      extras: undefined,
      error: 'Expected Rbracket, got: EOF',
      schema: ''
    }
  ]
  // Object projection is not supported, as we cannot tell ES which properties
  // to fetch. If the JMESPath query is made by the API, then this might work as
  // well.
  // [
  //   'object projection',
  //   'results.properties.electronic.*.type',
  //   {
  //     quantity: '?',
  //     path: '?',
  //     extras: [],
  //     error: undefined
  //   }
  // ],
])('parseJMESPath: %s ', (id, input, output) => {
  expect(isEqual(parseJMESPath(input), output)).toBe(true)
})
