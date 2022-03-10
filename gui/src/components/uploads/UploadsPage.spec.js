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
import 'regenerator-runtime/runtime'
import {
  render,
  screen,
  startAPI,
  closeAPI
} from '../conftest'
import UploadPage from './UploadPage'

test('Render upload page: unauthenticated', async () => {
  startAPI('tests.states.multiple_entries.multiple_entries', 'tests/data/multiple_entries/multiple_entries')
  render(<UploadPage uploadId={'dft_upload1'}/>)

  // Wait to load the page, i.e. wait for some text to appear
  await screen.findByText('unnamed upload')

  expect(screen.queryByTestId('step-prepare-and-upload-your-files')).toBeInTheDocument()
  expect(screen.queryByTestId('step-process-data')).toBeInTheDocument()
  expect(screen.queryByTestId('step-edit-author-metadata')).toBeNull()
  expect(screen.queryByTestId('step-publish')).toBeNull()

  closeAPI()
})
