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
import { screen, render } from '../conftest.spec'
import DeleteUploadsButton from './DeleteUploadsButton'
import userEvent from '@testing-library/user-event'

describe('test deleting different types of uploads', function() {
  test.each([
    [
      'single upload',
      'Delete selected upload',
      'Are you sure you want to delete the selected upload?',
      [{upload_id: 'a'}]
    ],
    [
      'multiple uploads',
      'Delete selected uploads',
      'You have selected 2 uploads. Are you sure you want to delete the selected uploads?',
      [{upload_id: 'a'}, {upload_id: 'b'}]
    ],
    [
      'upload with reviewer',
      'Delete selected upload',
      'The upload you are about to delete has been shared with at least one coauthor or reviewer. Are you sure you want to delete the selected upload?',
      [{upload_id: 'a', reviewers: [{user_id: 'b'}]}]
    ],
    [
      'upload with coauthor',
      'Delete selected upload',
      'The upload you are about to delete has been shared with at least one coauthor or reviewer. Are you sure you want to delete the selected upload?',
      [{upload_id: 'a', coauthors: [{user_id: 'b'}]}]
    ]
  ])('%s', async (name, tooltip, dialogText, uploads) => {
      const callback = jest.fn()
      render(<DeleteUploadsButton uploads={uploads} onConfirm={callback}/>)

      // Try clicking the button
      const button = screen.getByTooltip(tooltip)
      await userEvent.click(button)

      // See if text is correct
      expect(screen.getByText(dialogText)).toBeInTheDocument()

      // See that confirmation works
      const confirmButton = screen.getByText('Delete')
      expect(callback).not.toHaveBeenCalled()
      await userEvent.click(confirmButton)
      expect(callback).toHaveBeenCalled()
    }
  )
})
