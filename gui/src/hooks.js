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
import { useCallback, useEffect, useState, useMemo } from 'react'
import { debounce } from 'lodash'
import { useApi } from './components/api'

/**
 * Function for returning the current window size.
 *
 * @return {string} The number in approximated format.
 */
function getWindowSize() {
  return {windowWidth: window.innerWidth, windowHeight: window.innerHeight}
}

/**
 * Hook for returning the current window size.
 *
 * @return {string} The number in approximated format.
 */
export function useWindowSize() {
  const [size, setSize] = useState(getWindowSize())

  useEffect(() => {
    function handleResize() {
      setSize(getWindowSize())
    }

    window.addEventListener('resize', handleResize)
    return () => window.removeEventListener('resize', handleResize)
  }, [])

  return size
}

/**
 * Convenience hook for a simple boolean state.
 *
 * @param {bool} initialvalue The initial boolean value.
 * @return {array} Array containing current value and functions for setting to
 * true or false.
 */
export function useBoolState(initialValue) {
  const [value, setValue] = useState(initialValue)
  const setTrue = useCallback(() => setValue(true), [])
  const setFalse = useCallback(() => setValue(false), [])

  return [value, setTrue, setFalse]
}

/**
 * Hook for fetching suggestions for metainfo quantities. Notice that not all
 * quantities support suggestions.
 *
 * @param {array} quantities List of quantity names for which suggestions are fetched for.
 * @param {str} input Text input used for the suggestions. Must be non-empty in
 * order to fetch any suggestions.
 * @param {number} debounceTime The debounce time for restricting the number of
 * API calls.
 * @return {object} Array of suggestions contaiing the suggested value and the
 * quantity name.
 */
export function useSuggestions(quantities, input, debounceTime = 200, disabled = false) {
  const {api} = useApi()
  const [loading, setLoading] = useState(false)
  const [suggestions, setSuggestions] = useState([])

  // Used to retrieve suggestions for this field.
  const fetchSuggestions = useCallback((quantities, input) => {
    setLoading(true)
    api.suggestions(quantities, input)
      .then(data => {
        let res = []
        for (let quantity of quantities) {
          const esSuggestions = data[quantity]
          if (esSuggestions) {
            res = res.concat(esSuggestions.map(suggestion => ({
              category: quantity,
              value: suggestion.value
            })))
          }
        }
        setSuggestions(res)
      })
      .catch((error) => {
        console.log(error)
      })
      .finally(() => setLoading(false))
  }, [api])

  // Debounced version
  const fetchSuggestionsDebounced = useCallback(
    debounce(fetchSuggestions, debounceTime),
    [fetchSuggestions]
  )

  // Whenever input or the targeted quantities change, fetch suggestions
  useEffect(() => {
    const trimmedInput = input?.trim()
    if (trimmedInput && trimmedInput.length > 0) {
      fetchSuggestionsDebounced(quantities, input)
    } else {
      setSuggestions(old => old.length > 0 ? [] : old)
    }
  }, [fetchSuggestionsDebounced, quantities, input])

  const results = useMemo(() => [suggestions, loading], [suggestions, loading])
  return results
}
