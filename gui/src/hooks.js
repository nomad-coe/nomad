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
import { useCallback, useEffect, useState, useRef, useMemo } from 'react'
import { debounce, size } from 'lodash'
import { useApi } from './components/api'
import { staticSuggestions } from './components/search/FilterRegistry'

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
 * Hook for fetching suggestions for metainfo quantities. Uses a mixture of
 * metainfo Enumeration values and dynamic suggestions fetched through the API.
 * The fetching is debounced and retrieval of only the latest results is
 * ensured. Notice that not all quantities support suggestions.
 *
 * @param {array} quantities List of quantity names for which suggestions are
 *   fetched for. Notice that the suggestions are returned in this order.
 * @param {str} input Text input used for the suggestions. Must be non-empty in
 * order to fetch any suggestions.
 * @param {number} minLength Minimum input length before any dynamic
 *   suggestions are fetched through the API.
 * @param {number} suggestionsStatic An object containing any fixed suggestion
 *   values that should be added. Names that match a metainfo quantity will
 *   override any dynamical values.
 * @param {number} debounceTime The debounce time for restricting the number of
 * API calls.
 * @return {object} Array of suggestions containing the suggested value and the
 * quantity name.
 */
export function useSuggestions(quantities, input, minLength = 2, debounceTime = 150) {
  const {api} = useApi()
  const [loading, setLoading] = useState(false)
  const [suggestions, setSuggestions] = useState([])
  const currentRequest = useRef()

  // Used to retrieve suggestions for this field.
  const fetchSuggestions = useCallback((quantities, input) => {
    setLoading(true)

    // Used to flatten suggestions into a list that is correctly ordered.
    function flatten(value) {
      let suggestionsList = []
      quantities.forEach(q => {
        const options = value[q]
        if (size(options)) {
          suggestionsList = suggestionsList.concat(options)
        }
      })
      return suggestionsList
    }

    // Split the requested quantities into static and dynamic
    const quantitiesFixed = []
    const quantitiesDynamic = []
    for (const quantity of quantities) {
      if (quantity in staticSuggestions) {
        quantitiesFixed.push(quantity)
      } else {
        quantitiesDynamic.push(quantity)
      }
    }

    // Gather the fixed suggestions first
    const suggestionsTemp = {}
    for (const quantity of quantitiesFixed) {
      const fixed = staticSuggestions[quantity]
      suggestionsTemp[quantity] = fixed.filter(input)
    }

    // Start loading the dynamic suggestions
    if (input.length >= minLength) {
      const promise = api.suggestions(quantitiesDynamic, input)
      currentRequest.current = promise
      promise
        .then(data => {
          if (currentRequest.current === promise) {
            for (let quantity of quantitiesDynamic) {
              const esSuggestions = data[quantity]
              if (esSuggestions) {
                suggestionsTemp[quantity] = esSuggestions.map(suggestion => ({
                  category: quantity,
                  value: suggestion.value,
                  text: `${quantity}=${suggestion.value}`
                }))
              }
            }
            setSuggestions(flatten(suggestionsTemp))
          }
        })
        .catch((error) => {
          console.log(error)
        })
        .finally(() => {
          if (currentRequest.current === promise) {
            setLoading(false)
          }
        })
    } else {
      setSuggestions(flatten(suggestionsTemp))
      setLoading(false)
    }
  }, [api, minLength])

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
