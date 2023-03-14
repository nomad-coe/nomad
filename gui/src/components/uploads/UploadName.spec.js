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
import UploadName from './UploadName'
import userEvent from '@testing-library/user-event'
import { fireEvent } from '@testing-library/react'

describe('test name rendering', function() {
  test.each([
    ['undefined', undefined],
    ['null', null],
    ['empty string', ''],
    ['valid string', 'My Upload Name']
  ])('%s', async (name, uploadName) => {
      renderNoAPI(<UploadName upload_name={uploadName}/>)
      expect(screen.getByText(uploadName || 'unnamed upload')).toBeInTheDocument()
    }
  )
})

describe('test editing', function() {
  test.each([
    ['empty value', undefined, ''],
    ['empty string', '', ''],
    ['valid string', 'New Upload Name', 'New Upload Name'],
    ['test value trimming', ' New Name with Whitespaces  ', 'New Name with Whitespaces']
  ])('%s', async (name, input, expected) => {
      const onChange = jest.fn()
      renderNoAPI(<UploadName onChange={onChange}/>)

      // Press edit button
      const editButton = screen.getByRole('button')
      await userEvent.click(editButton)

      // Change value
      const inputField = screen.getByRole('textbox')
      await fireEvent.change(inputField, { target: { value: input } })

      // Press save button
      const saveButton = screen.getByRole('button')
      await userEvent.click(saveButton)

      // Expect the callback to be triggered with correct value
      expect(onChange).toHaveBeenCalledWith(expected)

      // Expect that the trimmed value is shown when re-entering edit mode
      await userEvent.click(screen.getByRole('button'))
      await screen.findByDisplayValue(expected, {exact: true, normalizer: (value) => value})
    }
  )
})
