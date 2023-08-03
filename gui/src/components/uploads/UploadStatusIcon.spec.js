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
import { screen, renderNoAPI } from '../conftest.spec'
import UploadStatusIcon from './UploadStatusIcon'

describe('test different states', function() {
  test.each([
    ['published', 'Published and accessible by everyone', {published: true}, {}],
    ['published, embargo, main author', 'Published with embargo by you and only accessible by you, coauthors and reviewers', {published: true, with_embargo: true, main_author: 'a'}, {sub: 'a'}],
    ['published, embargo, coauthor', 'Published with embargo and accessible by you as a coauthor', {published: true, with_embargo: true, main_author: 'a', coauthors: ['b']}, {sub: 'b'}],
    ['published, embargo, author', 'Published with embargo and accessible by you as a coauthor', {published: true, with_embargo: true, main_author: 'a', authors: [{user_id: 'b'}]}, {sub: 'b'}],
    ['published, embargo, reviewer', 'Published with embargo and accessible by you as a reviewer', {published: true, with_embargo: true, main_author: 'a', reviewers: ['b']}, {sub: 'b'}],
    ['published, embargo, viewer', 'Published with embargo and accessible by you as a reviewer', {published: true, with_embargo: true, main_author: 'a', viewers: [{user_id: 'b'}]}, {sub: 'b'}],
    ['published, embargo, external', 'Published with embargo and not accessible by you', {published: true, with_embargo: true, main_author: 'a', viewers: [{user_id: 'b'}]}, {sub: 'c'}],
    ['published, embargo, no user data', 'Published with embargo and might become accessible after login', {published: true, with_embargo: true, main_author: 'a', viewers: [{user_id: 'b'}]}, undefined],
    ['unpublished, main author', 'Unpublished, only accessible by you, coauthors and reviewers', {published: false, main_author: 'a'}, {sub: 'a'}],
    ['unpublished, coauthor', 'Unpublished, accessible by you as a coauthor', {published: false, main_author: 'a', coauthors: ['b']}, {sub: 'b'}],
    ['unpublished, author', 'Unpublished, accessible by you as a coauthor', {published: false, main_author: 'a', authors: [{user_id: 'b'}]}, {sub: 'b'}],
    ['unpublished, reviewer', 'Unpublished, accessible by you as a reviewer', {published: false, main_author: 'a', reviewers: ['b']}, {sub: 'b'}],
    ['unpublished, viewer', 'Unpublished, accessible by you as a reviewer', {published: false, main_author: 'a', viewers: [{user_id: 'b'}]}, {sub: 'b'}],
    ['unpublished, external', 'Unpublished', {published: false, main_author: 'a', viewers: [{user_id: 'b'}]}, {sub: 'c'}],
    ['unpublished, no user data', 'Unpublished', {published: false, main_author: 'a', viewers: [{user_id: 'b'}]}, undefined],
    ['no data', 'Upload status not available', undefined, undefined]
  ])('%s', async (name, tooltip, data, user) => {
      renderNoAPI(<UploadStatusIcon data={data} user={user}/>)
      expectUploadStatusIcon(tooltip)
    }
  )
})

/**
 * Asserts that an UploadStatusIcon is displayed correctly.
 *
 * @param {object} upload The upload object.
 * @param {string} state Expected state of the icon.
 */
export async function expectUploadStatusIcon(tooltip, root = screen) {
  expect(root.getByTooltip(tooltip)).toBeInTheDocument()
}
