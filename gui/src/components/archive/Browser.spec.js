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
import { join, basename } from 'path'
import React from 'react'
import PropTypes from 'prop-types'
import { waitFor } from '@testing-library/dom'
import { blockConsoleOutput, expectNoConsoleOutput, filteredConsoleOutput, consoleSpies, render, screen, within } from '../conftest.spec'
import { itemsInTreePath, checkLanes, navigateTo, selectItemAndWaitForRender, getLane } from './conftest.spec'
import Browser, { Item, Content, Compartment, Adaptor, Title, laneErrorBoundryMessage } from './Browser'

function checkTestLane({lane, laneIndex, lanePath, lastSegment, browserTree, rootTitle}) {
  expect(within(lane).getByText(laneIndex === 0 ? rootTitle : lastSegment)).toBeVisible() // Lane title

  itemsInTreePath(browserTree, lanePath).forEach(item => {
    expect(within(lane).getByText(item)).toBeVisible()
  })
}

const browserTree = {
  '': {cb: checkTestLane},
  'dir1': {cb: checkTestLane},
  'dir1/success': {cb: checkTestLane},
  'dir1/fail': {cb: checkTestLane}
}

class TestAdaptor extends Adaptor {
  constructor(path, title) {
    super()
    this.path = path
    this.title = title
  }

  itemAdaptor(key) {
    return new TestAdaptor(join(this.path, key))
  }

  render() {
    return <TestContent path={this.path} title={this.title} />
  }
}

function TestContent({path, title}) {
  if (path.endsWith('/fail')) {
    throw new Error('mocked render error')
  }
  return <Content key={path}>
    <Title title={title || basename(path)} label="test" />
    <Compartment>
      {
        itemsInTreePath(browserTree, path).map(itemKey => (
          <Item itemKey={itemKey} key={join(path, itemKey)}>
            {itemKey}
          </Item>
        ))
      }
    </Compartment>
  </Content>
}
TestContent.propTypes = {
  path: PropTypes.string.isRequired,
  title: PropTypes.string
}

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
    await expect(selectItemAndWaitForRender(getLane(1), 1, 'fail')).rejects.toThrow()
    expect(within(getLane(2)).queryByText(laneErrorBoundryMessage)).not.toBeNull()
    expect(filteredConsoleOutput().length).not.toBe(0)
  } finally {
    consoleSpies.logSpy.mockRestore()
    consoleSpies.errorSpy.mockRestore()
  }
})
