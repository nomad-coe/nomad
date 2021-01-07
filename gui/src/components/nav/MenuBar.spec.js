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
import { toHaveStyle, toBeVisible } from '@testing-library/jest-dom'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import { MenuBar, MenuBarItem, MenuBarMenu } from './MenuBar'
import { MemoryRouter } from 'react-router-dom'

function toBePrimaryColored(htmlElement) {
  return toHaveStyle(htmlElement, 'color: rgb(63, 81, 181)')
}

expect.extend({ toHaveStyle, toBeVisible, toBePrimaryColored })

describe('<NestedTopNav />', () => {
  const createExampleMenu = () => <MemoryRouter>
    <MenuBar selected="foo/foo-two">
      <MenuBarMenu name="foo">
        <MenuBarItem name="foo-one" route="/route" />
        <MenuBarItem name="foo-two" route="/route" />
      </MenuBarMenu>
      <MenuBarMenu name="bar" >
        <MenuBarItem name="bar-one" route="/route" />
        <MenuBarItem name="bar-two" route="/route" />
      </MenuBarMenu>
    </MenuBar>
    <div data-testid="content" />
  </MemoryRouter>

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

  it('opens menu on button click', async () => {
    render(createExampleMenu())

    expect(screen.getAllByRole('button').length).toBe(2)
    expect(screen.getByText('Foo-one')).not.toBeVisible()
    expect(screen.getByText('Bar-one')).not.toBeVisible()
    fireEvent.click(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )

    expect(screen.getByText('Foo-one')).toBeVisible()
    expect(screen.getByText('Bar-one')).not.toBeVisible()
  })

  it('changes open menu on hover', async () => {
    render(createExampleMenu())

    fireEvent.click(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )
    fireEvent.mouseOver(screen.getByText('Bar'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )

    expect(screen.getByText('Foo-one')).not.toBeVisible()
    expect(screen.getByText('Bar-one')).toBeVisible()
  })

  it('closes menu on item click', async () => {
    render(createExampleMenu())

    fireEvent.click(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )
    fireEvent.click(screen.getByText('Foo-one'))
    await waitFor(() =>
      expect(screen.getByText('Foo-one')).not.toBeVisible()
    )
  })

  it('closes menu on click somewhere', async () => {
    render(createExampleMenu())

    fireEvent.click(screen.getByText('Foo'))
    await waitFor(() =>
      expect(screen.getByRole('menu')).toBeVisible()
    )
    fireEvent.click(screen.getByTestId('content'))
    await waitFor(() =>
      expect(screen.getByText('Foo-one')).not.toBeVisible()
    )
  })
})
