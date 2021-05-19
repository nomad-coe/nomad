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
import { useCallback, useEffect } from 'react'
import { atom, useRecoilState } from 'recoil'

export const filtersState = atom({
  key: 'filters',
  default: new Set()
})

export function useFilters() {
  const [elements, setElements] = useRecoilState(filtersState)

  const addElement = useCallback(element => {
    setElements(old => {
      old.add(element)
      return old
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  const removeElement = useCallback(element => {
    setElements(old => {
      old.delete(element)
      return old
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  // Trigger search when query has changed. The querying is throttled so that
  // quick successive queries are bundled.
  useEffect(() => {
  }, [elements])

  return {
    elements: elements,
    setElements: setElements,
    addElement: addElement,
    removeElement: removeElement
  }
}
