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
import { fireEvent, waitFor } from '@testing-library/react'
import { render, screen } from '../conftest.spec'
import { MenuBar, MenuBarItem, MenuBarMenu } from './MenuBar'
import { Router } from 'react-router-dom'
import {createMemoryHistory} from 'history'

function toBePrimaryColored(received) {
  if (this.isNot) {
    expect(received).not.toEqual(expect.toHaveStyle('color: rgb(63, 81, 181)'))
  } else {
    expect(received).toEqual(expect.toHaveStyle('color: rgb(63, 81, 181)'))
  }
  return {pass: !this.isNot}
}

expect.extend({toBePrimaryColored})

describe('<NestedTopNav />', () => {
  const history = createMemoryHistory()
  history.push('/foo/foo-two')
  const createExampleMenu = () => <Router history={history}>
    <MenuBar>
      <MenuBarMenu label="Foo" route="/foo">
        <MenuBarItem label="Foo-one" route="/foo/foo-one" />
        <MenuBarItem label="Foo-two" route="/foo/foo-two" />
      </MenuBarMenu>
      <MenuBarMenu label="Bar" route="/bar" >
        <MenuBarItem label="Bar-one" route="/bar/bar-one" />
        <MenuBarItem label="Bar-two" route="/bar/bar-two" />
      </MenuBarMenu>
    </MenuBar>
    <div data-testid="content" />
  </Router>

  it('testid', () => {
    render(<div data-testid="myid">Hello</div>)
    expect(screen.getByTestId('myid')).toBeInTheDocument()
  })

  it('renders', () => {
    render(createExampleMenu())
    expect(screen.getAllByRole('button').length).toBe(2)
    expect(screen.getByText('Foo-one')).not.toBeVisible()
    expect(screen.getByText('Bar-one')).not.toBeVisible()
  })

  it('highlights selected', () => {
    render(createExampleMenu())
    expect(screen.getByText('Foo').parentElement).toBePrimaryColored()
    expect(screen.getByText('Foo-two').parentElement.parentElement).toBePrimaryColored()
    expect(screen.getByText('Foo-one').parentElement.parentElement).not.toBePrimaryColored()
    expect(screen.getByText('Bar').parentElement).not.toBePrimaryColored()
  })

  it('opens menu on mouse enter', async () => {
    render(createExampleMenu())
    expect(screen.getAllByRole('button').length).toBe(2)
    expect(screen.getByText('Foo-one')).not.toBeVisible()
    expect(screen.getByText('Bar-one')).not.toBeVisible()
    fireEvent.mouseEnter(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )

    expect(screen.getByText('Foo-one')).toBeVisible()
    expect(screen.getByText('Bar-one')).not.toBeVisible()
  })

  it('changes open menu when move from one to another', async () => {
    render(createExampleMenu())

    fireEvent.mouseEnter(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )
    fireEvent.mouseLeave(screen.getByText('Foo'))
    fireEvent.mouseEnter(screen.getByText('Bar'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )

    expect(screen.getByText('Foo-one')).not.toBeVisible()
    expect(screen.getByText('Bar-one')).toBeVisible()
  })

  it('closes menu on item click', async () => {
    render(createExampleMenu())

    fireEvent.mouseEnter(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )
    fireEvent.click(screen.getByText('Foo-one'))
    await waitFor(() =>
      expect(screen.getByText('Foo-one')).not.toBeVisible()
    )
  })

  it('closes menu on mouse leave', async () => {
    render(createExampleMenu())

    fireEvent.mouseEnter(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )
    fireEvent.mouseLeave(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByText('Foo-one')).not.toBeVisible()
    )
  })
})
