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
import { renderHook, act } from '@testing-library/react-hooks'
import { UnitProvider, useUnitContext } from './UnitContext'

const wrapper = ({ children }) => <UnitProvider
  initialUnitSystems={{
    SI: {units: {'length': {'definition': 'meter'}}},
    AU: {units: {'length': {'definition': 'ang'}}}
  }}
  initialSelected='SI'
>
  {children}
</UnitProvider>

test('the initial selection is returned correctly', () => {
  const { result } = renderHook(() => useUnitContext(), { wrapper })
  expect(result.current.selected).toBe('SI')
})

test('updating unit system selection works', () => {
  const { result } = renderHook(() => useUnitContext(), { wrapper })
  act(() => {
    result.current.setSelected('AU')
  })
  expect(result.current.selected).toBe('AU')
})

test('the initial units are returned correctly', () => {
  const { result } = renderHook(() => useUnitContext(), { wrapper })
  expect(result.current.units.length.definition).toBe('meter')
})

test('updating units works', () => {
  const { result } = renderHook(() => useUnitContext(), { wrapper })
  act(() => {
    result.current.setUnits(old => ({...old, length: {definition: 'mm'}}))
  })
  expect(result.current.units.length.definition).toBe('mm')
})

test('switching between scopes', () => {
  const { result, rerender } = renderHook(
    ({ scopeKey, defaultUnitSystem }) => useUnitContext(scopeKey, defaultUnitSystem),
    {
      wrapper: wrapper,
      initialProps: {
        scopeKey: 'schemaId',
        defaultUnitSystem: 'AU'
      }
    }
  )
  // In first render of a schema the default unit system should appear
  expect(result.current.selected).toBe('AU')
  expect(result.current.units.length.definition).toBe('ang')
  // In schema the scope `Schema` should be selected as default
  expect(result.current.scope).toBe('Schema')

  // Go to another scope i.e. search
  rerender({
    scopeKey: undefined,
    defaultUnitSystem: undefined
  })
  expect(result.current.selected).toBe('SI')
  expect(result.current.units.length.definition).toBe('meter')
  expect(result.current.scope).toBe(undefined)

  // Change in schema scope
  rerender({
    scopeKey: 'schemaId',
    defaultUnitSystem: 'AU'
  })
  act(() => {
    result.current.setSelected('SI')
  })
  expect(result.current.selected).toBe('SI')
  expect(result.current.units.length.definition).toBe('meter')
  act(() => {
    result.current.setUnits(old => ({...old, length: {definition: 'mm'}}))
  })
  expect(result.current.units.length.definition).toBe('mm')

  // Change in global scope
  rerender({
    scopeKey: undefined,
    defaultUnitSystem: undefined
  })
  expect(result.current.selected).toBe('SI')
  act(() => {
    result.current.setSelected('AU')
  })
  expect(result.current.selected).toBe('AU')
  expect(result.current.units.length.definition).toBe('ang')
  expect(result.current.scope).toBe(undefined)

  // Test the last units in schema scope
  rerender({
    scopeKey: 'schemaId',
    defaultUnitSystem: 'AU'
  })
  expect(result.current.selected).toBe('SI')
  expect(result.current.units.length.definition).toBe('mm')
  expect(result.current.scope).toBe('Schema')
  act(() => {
    result.current.setScope('Global')
  })
  expect(result.current.scope).toBe('Global')
  expect(result.current.selected).toBe('AU')
  expect(result.current.units.length.definition).toBe('ang')

  // Test the last units in global scope
  rerender({
    scopeKey: undefined,
    defaultUnitSystem: undefined
  })
  expect(result.current.selected).toBe('AU')
  expect(result.current.units.length.definition).toBe('ang')
})
