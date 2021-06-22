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
import {
  atom,
  atomFamily,
  selector,
  useSetRecoilState,
  useRecoilValue,
  useRecoilState,
  useRecoilCallback
} from 'recoil'
import _ from 'lodash'
import { useApi } from '../apiV1'
import { setToArray } from '../../utils'
import { Quantity } from '../../units'

/**
 * Each search quantity is here mapped into a separate Recoil.js Atom. This
 * allows components to hook into individual search parameters (both for setting
 * and reading their value). This performs much better than having one large
 * Atom for the entire query, as this would cause all of the hooked components
 * to render even if they are not affected by some other search quantity.
 * Re-renders became problematic with large and complex components (e.g. the
 * periodic table), for which the re-rendering takes significant time. Another
 * approach would have been to try and Memoize each sufficiently complex
 * component, but this quickly becomes a hard manual task.
 */

export const filtersElements = [
  'results.material.elements',
  'results.material.chemical_formula_hill',
  'results.material.chemical_formula_anonymous',
  'results.material.n_elements'
]

export const filtersMaterial = [
  'results.material.structural_type',
  'results.material.functional_type',
  'results.material.compound_type',
  'results.material.material_name'
]

export const filtersElectronic = [
  'results.properties.electronic.band_structure_electronic.channel_info.band_gap',
  'results.properties.electronic.band_structure_electronic.channel_info.band_gap_type'
]

export const filtersSymmetry = [
  'results.material.symmetry.bravais_lattice',
  'results.material.symmetry.crystal_system',
  'results.material.symmetry.hall_symbol',
  'results.material.symmetry.point_group',
  'results.material.symmetry.space_group_symbol',
  'results.material.symmetry.prototype_aflow_id',
  'results.material.symmetry.structure_name',
  'results.material.symmetry.strukturbericht_designation'
]

export const filtersMethod = [
  'results.method.method_name',
  'results.method.simulation.program_name',
  'results.method.simulation.program_version'
]

export const filtersDFT = [
  'results.method.simulation.dft.basis_set_type',
  'results.method.simulation.dft.basis_set_name',
  'results.method.simulation.dft.core_electron_treatment',
  'results.method.simulation.dft.van_der_Waals_method',
  'results.method.simulation.dft.relativity_method',
  'results.method.simulation.dft.smearing_type'
]

export const filtersGW = [
  'results.method.simulation.gw.gw_type'
]

export const filtersAuthor = [
  'authors.name',
  'external_db'
]

export const filtersDataset = [
  'datasets.name',
  'datasets.doi'
]

export const filtersIDs = [
  'results.material.material_id',
  'upload_id',
  'entry_id',
  'datasets.dataset_id'
]

export let filtersAll = []
filtersAll = filtersAll.concat(filtersElements)
filtersAll = filtersAll.concat(filtersMaterial)
filtersAll = filtersAll.concat(filtersElectronic)
filtersAll = filtersAll.concat(filtersSymmetry)
filtersAll = filtersAll.concat(filtersMethod)
filtersAll = filtersAll.concat(filtersDFT)
filtersAll = filtersAll.concat(filtersGW)
filtersAll = filtersAll.concat(filtersAuthor)
filtersAll = filtersAll.concat(filtersDataset)
filtersAll = filtersAll.concat(filtersIDs)

export const queryFamily = atomFamily({
  key: 'queryFamily',
  default: undefined
})

let index = 0

// Menu open state
export const menuOpen = atom({
  key: 'isMenuOpen',
  default: false
})
export function useMenuOpenState() {
  return useRecoilState(menuOpen)
}
export function useSetMenuOpen() {
  return useSetRecoilState(menuOpen)
}

// Exclusive search state
export const exclusive = atom({
  key: 'exclusive',
  default: false
})
export function useExclusive() {
  return useRecoilValue(exclusive)
}
export function useExclusiveState() {
  return useRecoilState(exclusive)
}
export function useSetExclusive() {
  return useSetRecoilState(exclusive)
}

/**
 * Returns a function that can be called to reset all current filters.
 */
export function useResetFilters() {
  const reset = useRecoilCallback(({reset}) => () => {
    for (let filter of filtersAll) {
      reset(queryFamily(filter))
    }
  }, [])
  return reset
}

/**
 * This hook will expose a function for reading filter values for a specific
 * quantity. Use this hook if you intend to only view the filter values and are
 * not interested in setting the filter.
 *
 * @param {*} quantity Name of the quantity. Should exist in searchQuantities.json.
 * @returns currently set filter value.
 */
export function useFilterValue(quantity) {
  return useRecoilValue(queryFamily(quantity))
}
/**
 * This hook will expose a function for setting filter values for a specific
 * quantity. Use this hook if you intend to only set the filter value and are
 * not interested in the query results.
 *
 * @param {*} quantity Name of the quantity to set. Should exist in searchQuantities.json.
 * @param {Set} set An optional Set that keeps track of hooked filters.
 * @returns function for setting the value for the given quantity
 */
export function useSetFilter(quantity) {
  return useSetRecoilState(queryFamily(quantity))
}
/**
 * This hook will expose a function for getting and setting filter values for a
 * specific quantity. Use this hook if you intend to both read and write the
 * filter value.
 *
 * @param {*} quantity Name of the quantity to set. Should exist in searchQuantities.json.
 * @returns array containing the filter value and setter function for it.
 */
export function useFilterState(quantity) {
  return useRecoilState(queryFamily(quantity))
}

/**
 * This hook will expose a function for getting and setting filter values for
 * the specified list of quantities. Use this hook if you intend to both read
 * and write the filter values.
 *
 * @param {*} quantities Names of the quantities. Should exist in searchQuantities.json.
 * @param {string} id Unique ID for this set of Filters (needed by Recoil.js)
 * @returns array containing the filter value and setter function for it.
 */
export function useFiltersState(quantities) {
  // We dynamically create a Recoil.js selector that is subscribed to the
  // filters specified in the input. This way only the specified filters will
  // cause a render.

  // Recoil.js requires that each selector/atom has an unique id. Because this
  // hook can be called dynamically, we simply generate the ID sequentially.
  const id = `dynamic_selector${index}`
  index += 1
  const filterState = useMemo(() => {
    return selector({
      key: id,
      get: ({get}) => {
        const query = {}
        for (let key of quantities) {
          const filter = get(queryFamily(key))
          query[key] = filter
        }
        return query
      },
      set: ({set}, [key, value]) => {
        set(queryFamily(key), value)
      }
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return useRecoilState(filterState)
}

/**
 * This selector aggregates all the currently set filters into a single query
 * object used by the API.
 */
const queryState = selector({
  key: 'query',
  get: ({get}) => {
    let query = {}
    for (let key of filtersAll) {
      const filter = get(queryFamily(key))
      if (filter !== undefined) {
        query[key] = filter
      }
    }
    return query
  }
})
export function useQueryValue() {
  return useRecoilValue(queryState)
}

/**
 * Hook for returning the current search object.
 *
 * @returns {object} Object containing the search object.
 */
export function useSearch() {
  const query = useRecoilValue(queryState)
  const result = useMemo(() => {
    return {
      owner: 'all',
      query: query,
      pagination: {
        page: 1,
        page_size: 10,
        order: 'desc',
        order_by: 'upload_time'
      }
    }
  }, [query])
  return result
}

/**
 * Hook for retrieving the most up-to-date aggreagtion results for a specific
 * quantity, taking into account the current search context.
 *
 * @param {string} quantity
 * @param {string} type
 * @param {bool} restrict
 * @param {bool} update Whether the hook needs to react to changes in the
 * current query context. E.g. if the component showing the data is not visible,
 * this can be set to false.
 * @returns {array} The data-array returned by the API.
 */
export function useAgg(quantity, type, restrict = false, update = true, delay = 200) {
  const api = useApi()
  const [results, setResults] = useState()
  const query = useQueryValue()
  const exclusive = useExclusive()
  const firstRender = useRef(true)

  // Pretty much all of the required pre-processing etc. should be done in this
  // function, as it is the final one that gets called after the debounce
  // interval.
  const apiCall = useCallback((query, exclusive) => {
    // If the restrict option is enabled, the filters targeting the specified
    // quantity will be removed. This way all possible options pre-selection can
    // be returned.
    let queryCleaned = {...query}
    if (restrict && query && quantity in query) {
      queryCleaned[quantity] = undefined
    }
    queryCleaned = cleanQuery(queryCleaned, exclusive)

    const aggs = {
      [quantity]: {
        [type]: {
          quantity: quantity,
          size: 50
        }
      }
    }
    const search = {
      owner: 'visible',
      query: queryCleaned,
      aggregations: aggs,
      pagination: {page_size: 0},
      required: { include: [] }
    }

    api.queryEntry(search, false)
      .then(data => {
        setResults(data)
      })
  }, [api, quantity, restrict, type])

  // This is a debounced version of apiCall.
  const debounced = useCallback(_.debounce(apiCall, delay), [])

  // The API call is made immediately on first render. On subsequent renders it
  // will be debounced.
  useEffect(() => {
    if (!update && !firstRender.current) {
      return
    }
    if (firstRender.current) {
      apiCall(query, exclusive)
      firstRender.current = false
    } else {
      debounced(query, exclusive)
    }
  }, [apiCall, debounced, query, exclusive, update])

  return results && results.aggregations[quantity][type].data
}

/**
 * Hook for returning the current query results.
 *
 * @param {number} delay The debounce delay in milliseconds.
 * @returns {object} Object containing the search results under 'results' and
 * the used query under 'search'.
 */
export function useResults(delay = 400) {
  const api = useApi()
  const firstRender = useRef(true)
  const search = useSearch()
  const exclusive = useExclusive()
  const [results, setResults] = useState()

  // The results are fetched as a side effect in order to not block the
  // rendering. This causes two renders: first one without the data, the second
  // one with the data.
  const apiCall = useCallback((search, exclusive) => {
    const finalSearch = {...search}
    finalSearch.query = cleanQuery(finalSearch.query, exclusive)

    api.queryEntry(finalSearch)
      .then(data => {
        setResults(data)
      })
  }, [api])

  // This is a debounced version of apiCall.
  const debounced = useCallback(_.debounce(apiCall, delay), [])

  // The API call is made immediately on first render. On subsequent renders it
  // will be debounced.
  useEffect(() => {
    if (firstRender.current) {
      apiCall(search, exclusive)
      firstRender.current = false
    } else {
      debounced(search, exclusive)
    }
  }, [apiCall, debounced, search, exclusive])

  return results
}

// Converts all sets to arrays and convert all Quantities into their SI unit
// values
function cleanQuery(obj, exclusive) {
  let newObj = {}
  for (let [k, v] of Object.entries(obj)) {
    let newValue

    // Special handling for elements. There are two search modes: exclusive and
    // inclusive.
    if (k === 'results.material.elements') {
      if (v.size === 0) {
        continue
      }
      if (exclusive) {
        k = `${k}.exclusive`
        newValue = setToArray(v).sort().join(' ')
      } else {
        k = `${k}:any`
        newValue = setToArray(v)
      }
    } else {
      if (v instanceof Set) {
        newValue = setToArray(v)
        k = `${k}:any`
      } else if (v instanceof Quantity) {
        newValue = v.toSI()
      } else if (Array.isArray(v)) {
        newValue = v
        k = `${k}:any`
      } else if (typeof v === 'object' && v !== null) {
        newValue = cleanQuery(v, exclusive)
      } else {
        newValue = v
      }
    }
    newObj[k] = newValue
  }
  return newObj
}
