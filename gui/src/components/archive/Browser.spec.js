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
import { waitFor } from '@testing-library/dom'
import { blockConsoleOutput, expectNoConsoleOutput, filteredConsoleOutput, consoleSpies, render, screen, within } from '../conftest.spec'
import { checkLanes, navigateTo, selectItemAndWaitForRender, getLane, TestAdaptor, browserTree } from './conftest.spec'
import Browser, { laneErrorBoundryMessage } from './Browser'

test('Test browser lane error boundry', async () => {
  blockConsoleOutput()

  try {
    const browserConfig = { browserTree, rootTitle: 'Root title' }
    const adaptor = new TestAdaptor('', browserConfig.rootTitle)
    render(<Browser adaptor={adaptor} />)

    await waitFor(() => {
      expect(screen.getByText(browserConfig.rootTitle)).toBeVisible()
    })
    await checkLanes('', browserConfig)
    await navigateTo('dir1', browserConfig)
    await navigateTo('dir1/success', browserConfig)
    expectNoConsoleOutput()
    // Call lane which fails to render
    await expect(selectItemAndWaitForRender(1, 'fail')).rejects.toThrow()
    expect(within(getLane(2)).queryByText(laneErrorBoundryMessage)).not.toBeNull()
    expect(filteredConsoleOutput().length).not.toBe(0)
  } finally {
    consoleSpies.logSpy.mockRestore()
    consoleSpies.errorSpy.mockRestore()
  }
})
