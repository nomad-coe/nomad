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
import {
  debounce,
  isEmpty,
  isArray,
  isPlainObject,
  isNil
} from 'lodash'
import qs from 'qs'
import { useHistory } from 'react-router-dom'
import { useApi } from '../apiV1'
import { setToArray, formatMeta, parseMeta } from '../../utils'
import { Quantity } from '../../units'

let index = 0
export const quantities = new Set()
export const quantityGroups = new Map()
export const quantityAggregations = new Map()
export const quantityAbbreviations = new Map()
export const quantityFullnames = new Map()
export const labelMaterial = 'Material'
export const labelElements = 'Elements / Formula'
export const labelSymmetry = 'Symmetry'
export const labelMethod = 'Method'
export const labelDFT = 'DFT'
export const labelElectronic = 'Electronic'
export const labelAuthor = 'Author / Origin'
export const labelDataset = 'Dataset'
export const labelIDs = 'IDs'

/**
 * This function is used to register a quantity within the FilterContext.
 *
 * Only registered quantities may be searched for. The registration must happen
 * before any components use the filters associated with quantities. This is
 * because:
 *  - The initial aggregation results must be fetched before any components
 *  using the filter values are rendered.
 *  - Several components need to know the list of available filter (e.g. the
 *  search bar and  the search panel). If quantities are only registered during
 *  component initialization, it may already be too late to update other
 *  components.
 *
 * @param {string} quantity Name of the quantity. Should exist in searchQuantities.json.
 * @param {string} group The group into which the quantity belongs to. Groups
 * are used to e.g. in showing FilterSummaries about a group of filters.
 * @param {string} agg Possible aggregation associated with the quantity. If
 * provided, the initial aggregation value will be prefetched for this quantity.
 */
function registerQuantity(name, group, agg) {
  quantities.add(name)
  quantityGroups.has(group) ? quantityGroups.get(group).add(name) : quantityGroups.set(group, new Set([name]))
  if (agg) {
    quantityAggregations[name] = agg
  }
  const abbreviation = name.split('.').pop()
  const oldName = quantityAbbreviations.get(abbreviation)
  if (oldName === undefined) {
    quantityAbbreviations.set(name, abbreviation)
    quantityFullnames.set(abbreviation, name)
  } else {
    quantityFullnames.delete(abbreviation)
    quantityAbbreviations.set(name, name)
    quantityAbbreviations.set(oldName, oldName)
  }
}

registerQuantity('results.material.structural_type', labelMaterial, 'terms')
registerQuantity('results.material.functional_type', labelMaterial, 'terms')
registerQuantity('results.material.compound_type', labelMaterial, 'terms')
registerQuantity('results.material.material_name', labelMaterial)
registerQuantity('results.material.elements', labelElements, 'terms')
registerQuantity('results.material.elements_exclusive', labelElements, 'terms')
registerQuantity('results.material.chemical_formula_hill', labelElements)
registerQuantity('results.material.chemical_formula_anonymous', labelElements)
registerQuantity('results.material.n_elements', labelElements, 'min_max')
registerQuantity('results.material.symmetry.bravais_lattice', labelSymmetry, 'terms')
registerQuantity('results.material.symmetry.crystal_system', labelSymmetry, 'terms')
registerQuantity('results.material.symmetry.structure_name', labelSymmetry, 'terms')
registerQuantity('results.material.symmetry.strukturbericht_designation', labelSymmetry, 'terms')
registerQuantity('results.material.symmetry.space_group_symbol', labelSymmetry)
registerQuantity('results.material.symmetry.point_group', labelSymmetry)
registerQuantity('results.material.symmetry.hall_symbol', labelSymmetry)
registerQuantity('results.material.symmetry.prototype_aflow_id', labelSymmetry)
registerQuantity('results.method.method_name', labelMethod, 'terms')
registerQuantity('results.method.simulation.program_name', labelMethod, 'terms')
registerQuantity('results.method.simulation.program_version', labelMethod)
registerQuantity('results.method.simulation.dft.basis_set_type', labelDFT, 'terms')
registerQuantity('results.method.simulation.dft.core_electron_treatment', labelDFT, 'terms')
registerQuantity('results.method.simulation.dft.relativity_method', labelDFT, 'terms')
registerQuantity('results.properties.electronic.band_structure_electronic.channel_info.band_gap', labelElectronic, 'min_max')
registerQuantity('results.properties.electronic.band_structure_electronic.channel_info.band_gap_type', labelElectronic, 'terms')
registerQuantity('results.properties.available_properties', labelElectronic, 'terms')
registerQuantity('external_db', labelAuthor, 'terms')
registerQuantity('authors.name', labelAuthor)
registerQuantity('upload_time', labelAuthor, 'min_max')
registerQuantity('datasets.name', labelDataset)
registerQuantity('datasets.doi', labelDataset)
registerQuantity('entry_id', labelIDs)
registerQuantity('upload_id', labelIDs)
registerQuantity('results.material.material_id', labelIDs)
registerQuantity('datasets.dataset_id', labelIDs)

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
export const queryFamily = atomFamily({
  key: 'queryFamily',
  default: undefined
})

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

// Whether the search is initialized.
export const initialized = atom({
  key: 'initialized',
  default: false
})

// Exclusive search state
export const exclusive = atom({
  key: 'exclusive',
  default: false
})
export function useExclusive() {
  return useRecoilValue(exclusive)
}
export function useExclusiveState() {
  const value = useRecoilValue(exclusive)
  const setter = useSetExclusive()
  return [value, setter]
}
export function useSetExclusive() {
  const setter = useSetRecoilState(exclusive)
  const [elements, setElements] = useRecoilState(queryFamily('results.material.elements'))
  const [elementsEx, setElementsEx] = useRecoilState(queryFamily('results.material.elements_exclusive'))

  const set = useCallback((exclusive) => {
    setter(exclusive)
    if (exclusive) {
      setElements(new Set())
      setElementsEx(elements)
    } else {
      setElements(elementsEx)
      setElementsEx(new Set())
    }
  }, [elements, elementsEx, setElements, setElementsEx, setter])

  return set
}

/**
 * Returns a function that can be called to reset all current filters.
 *
 * @returns Function for resetting all filters.
 */
export function useResetFilters() {
  const reset = useRecoilCallback(({reset}) => () => {
    for (let filter of quantities) {
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
 * This Recoil.js selector aggregates all the currently set filters into a
 * single query object used by the API.
 */
const queryState = selector({
  key: 'query',
  get: ({get}) => {
    const inited = get(initialized)
    if (!inited) {
      return undefined
    }
    let query = {}
    for (let key of quantities) {
      const filter = get(queryFamily(key))
      if (filter !== undefined) {
        query[key] = filter
      }
    }
    return query
  },
  set: ({ get, set, reset }, data) => {
    set(initialized, true)
    for (let filter of quantities) {
      reset(queryFamily(filter))
    }
    if (data) {
      for (const [key, value] of Object.entries(data)) {
        set(queryFamily(key), value)
      }
    }
  }
})

/**
 * Hook used to initialize the query. Must be called once before any results can
 * be fetched. By default uses the URL query string to initialize the query.
*/
export function useInitQuery() {
  const setQuery = useSetRecoilState(queryState)
  const setInitialized = useSetRecoilState(initialized)

  // If a query string is available, initialize the query from it. No results
  // are returned by this hook until the setQuery has made its round back. This
  // prevents double renders.
  useEffect(() => {
    const location = window.location.href
    const queryString = location.split('?')
    const qs = queryString.length !== 1 && qsToQuery(queryString.pop())
    setInitialized(true)
    setQuery(qs || {})
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])
}

export function useQuery() {
  return useRecoilValue(queryState)
}

/**
 * Hook for writing a query object to the query string.
 *
 * @returns {object} Object containing the search object.
 */
export function useUpdateQueryString() {
  const history = useHistory()

  const updateQueryString = useCallback((query) => {
    const queryString = queryToQs(query)
    history.replace(history.location.pathname + '?' + queryString)
  }, [history])

  return updateQueryString
}

/**
 * Converts a query string into a valid query object.
 *
 * @param {string} queryString URL querystring, encoded or not.
 * @returns Returns an object containing the filters. Values are converted into
 * datatypes that are directly compatible with the filter components.
 */
function qsToQuery(queryString) {
  const query = qs.parse(queryString, { comma: true })
  const newQuery = {}
  for (let [key, value] of Object.entries(query)) {
    const split = key.split(':')
    key = split[0]
    let newKey = quantityFullnames.get(key) || key
    const {type, parser} = parseMeta(newKey)
    if (split.length !== 1) {
      const op = split[1]
      const oldValue = newQuery[newKey]
      if (!oldValue) {
        newQuery[newKey] = {[op]: parser(value)}
      } else {
        newQuery[newKey][op] = parser(value)
      }
    } else {
      if (isArray(value)) {
        value = new Set(value)
      } else if (isPlainObject(value)) {
        if (!isNil(value.gte)) {
          value.gte = parser(value.gte)
        }
        if (!isNil(value.lte)) {
          value.lte = parser(value.lte)
        }
      } else {
        value = parser(value)
        if (type !== 'number' && type !== 'timestamp') {
          value = new Set([value])
        }
      }
      newQuery[newKey] = value
    }
  }
  return newQuery
}

/**
 * Converts a query into a valid query string.
 * @param {object} query A query object representing the currently active
 * filters.
 * @returns URL querystring, not encoded if possible to improve readability.
 */
function queryToQs(query) {
  const newQuery = {}
  for (const [key, value] of Object.entries(query)) {
    const {formatter} = formatMeta(key, false)
    let newValue
    const newKey = quantityAbbreviations.get(key)
    if (isPlainObject(value)) {
      if (!isNil(value.gte)) {
        newQuery[`${newKey}:gte`] = formatter(value.gte)
      }
      if (!isNil(value.lte)) {
        newQuery[`${newKey}:lte`] = formatter(value.lte)
      }
    } else {
      if (isArray(value)) {
        newValue = value.map(formatter)
      } else if (value instanceof Set) {
        newValue = [...value].map(formatter)
      } else {
        newValue = formatter(value)
      }
      newQuery[newKey] = newValue
    }
  }
  return qs.stringify(newQuery, {indices: false, encode: false})
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

export const initialAggs = atom({
  key: 'initialAggs',
  default: undefined
})

/**
 * Hook for returning the current search object.
 *
 * @returns {object} Object containing the search object.
 */
export function useInitialAgg(quantity, agg) {
  const aggs = useRecoilValue(initialAggs)
  let data = aggs && aggs?.aggregations?.[quantity]?.[agg].data
  if (agg === 'min_max' && !data) {
    data = [undefined, undefined]
  }
  return data
}

export function useInitialAggs() {
  // Request initial aggregation values for all filter components that are on
  // the search page. This is only done once.
  const api = useApi()
  const setInitialAggs = useSetRecoilState(initialAggs)

  useEffect(() => {
    const aggRequest = {}
    for (const [quantity, agg] of Object.entries(quantityAggregations)) {
      aggRequest[quantity] = {
        [agg]: {
          quantity: quantity,
          size: 500
        }
      }
    }
    // Make the request
    const search = {
      owner: 'all',
      query: {},
      aggregations: aggRequest,
      pagination: {page_size: 0}
    }

    api.queryEntry(search, false)
      .then(data => {
        setInitialAggs(data)
      })
  }, [api, setInitialAggs])
}

/**
 * Hook for retrieving the most up-to-date aggregation results for a specific
 * quantity, taking into account the current search context.
 *
 * @param {string} quantity The quantity name
 * @param {string} type ElasticSearch aggregation type
 * @param {bool} restrict If true, the query filters targeting this particular
 * quantity will be removed. This makes it possible to return all possible
 * values for dropdowns etc.
 * @param {bool} update Whether the hook needs to react to changes in the
 * current query context. E.g. if the component showing the data is not visible,
 * this can be set to false.
 *
 * @returns {array} The data-array returned by the API.
 */
export function useAgg(quantity, type, restrict = false, update = true, delay = 200) {
  const api = useApi()
  const [results, setResults] = useState(type === 'min_max' ? [undefined, undefined] : undefined)
  const initialAgg = useInitialAgg(quantity, type)
  const query = useQuery()
  const exclusive = useExclusive()
  const firstLoad = useRef(true)

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
          size: 500
        }
      }
    }
    const search = {
      owner: 'all',
      query: queryCleaned,
      aggregations: aggs,
      pagination: {page_size: 0},
      required: { include: [] }
    }

    api.queryEntry(search, false)
      .then(data => {
        let cleaned = data && data.aggregations[quantity][type].data
        if (type === 'min_max' && !cleaned) {
          cleaned = [undefined, undefined]
        }
        firstLoad.current = false
        setResults(cleaned)
      })
  }, [api, quantity, restrict, type])

  // This is a debounced version of apiCall.
  const debounced = useCallback(debounce(apiCall, delay), [])

  // The API call is made immediately on first render. On subsequent renders it
  // will be debounced.
  useEffect(() => {
    if (!update) {
      return
    }
    if (firstLoad.current) {
      // Fetch the initial aggregation values if no query
      // is specified.
      if (isEmpty(query)) {
        setResults(initialAgg)
      // Make an immediate request for the aggregation values if query has been
      // specified.
      } else {
        apiCall(query, exclusive)
      }
    } else {
      debounced(query, exclusive)
    }
  }, [apiCall, debounced, query, exclusive, update, initialAgg])

  return results
}

/**
 * Hook for returning a set of results based on the currently set query together
 * with a function for retrieving a new set of results.
 *
 * @param {object} pagination The pagination settings as used by the API. Notice
 * that 'page', and 'page_after_value' will be ignored, as they are controlled
 * by the hook.
 * @param {bool} exclusive Whether to use exclusive element search.
 * @param {number} delay The debounce delay in milliseconds.
 *
 * @returns {object} Object containing the search results under 'results' and
 * the used query under 'search'.
 */
export function useScrollResults(pageSize, orderBy, order, exclusive, delay = 400) {
  const api = useApi()
  const firstRender = useRef(true)
  const [results, setResults] = useState()
  const pageNumber = useRef(1)
  const query = useQuery(true)
  const updateQueryString = useUpdateQueryString()
  const pageAfterValue = useRef()
  const searchRef = useRef()
  const loading = useRef(false)
  const total = useRef(0)

  // The results are fetched as a side effect in order to not block the
  // rendering. This causes two renders: first one without the data, the second
  // one with the data.
  const apiCall = useCallback((query, pageSize, orderBy, order, exclusive) => {
    pageAfterValue.current = undefined
    const cleanedQuery = cleanQuery(query, exclusive)
    const search = {
      owner: 'all',
      query: cleanedQuery,
      pagination: {
        page_size: pageSize,
        order_by: orderBy,
        order: order,
        page_after_value: pageAfterValue.current
      }
    }
    searchRef.current = search

    loading.current = true
    api.queryEntry(search)
      .then(data => {
        pageAfterValue.current = data.pagination.next_page_after_value
        total.current = data.pagination.total
        setResults(data)
        loading.current = false
      })
  }, [api])

  // This is a debounced version of apiCall.
  const debounced = useCallback(debounce(apiCall, delay), [])

  // Used to load the next bath of results
  const next = useCallback(() => {
    if (loading.current) {
      return
    }
    pageNumber.current += 1
    searchRef.current.pagination.page_after_value = pageAfterValue.current
    loading.current = true
    api.queryEntry(searchRef.current)
      .then(data => {
        pageAfterValue.current = data.pagination.next_page_after_value
        total.current = data.pagination.total
        setResults(old => {
          data.data = old.data.concat(data.data)
          return data
        })
        loading.current = false
      })
  }, [api])

  // Whenever the query changes, we make a new query that resets pagination and
  // shows the first batch of results.
  useEffect(() => {
    // If the initial query is not yet ready, do nothing
    if (query === undefined) {
      return
    }
    if (firstRender.current) {
      apiCall(query, pageSize, orderBy, order, exclusive)
      firstRender.current = false
    } else {
      updateQueryString(query)
      debounced(query, pageSize, orderBy, order, exclusive)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [apiCall, debounced, query, exclusive, pageSize, order, orderBy])

  // Whenever the ordering changes, we perform a single API call that fetches
  // results in the new order. The amount of fetched results is based on the
  // already loaded amount.
  // TODO

  return {
    results: results,
    next: next,
    page: pageNumber.current,
    total: total.current
  }
}

/**
 * Converts all sets to arrays and convert all Quantities into their SI unit
 * values.
 *
 * Should only be called when making the final API call, as during the
 * construction of the query it is much more convenient to store filters within
 * e.g. Sets.
 *
 * @param {number} query The query object
 * @param {bool} exclusive The chemical element search mode.
 *
 * @returns {object} A copy of the object with certain items cleaned into a
 * format that is supported by the API.
 */
export function cleanQuery(query, exclusive) {
  let newQuery = {}
  for (let [k, v] of Object.entries(query)) {
    let newValue

    // If a regular elements query is made, we add the ':all'-prefix.
    if (k === 'results.material.elements') {
      if (v.size === 0) {
        continue
      }
      k = `${k}:all`
      newValue = setToArray(v)
    // If an exlusive elements query is made, we sort the elements and
    // concatenate them into a single string. This value we can then use to
    // target the special field reserved for exclusive queries. TODO: Maybe the
    // API could directly support an ':only'-postfix?
    } else if (k === 'results.material.elements_exclusive') {
      if (v.size === 0) {
        continue
      }
      newValue = setToArray(v).sort().join(' ')
    } else {
      if (v instanceof Set) {
        newValue = setToArray(v)
        if (newValue.length === 0) {
          newValue = undefined
        } else {
          newValue = newValue.map((item) => item instanceof Quantity ? item.toSI() : item)
        }
        k = `${k}:any`
      } else if (v instanceof Quantity) {
        newValue = v.toSI()
      } else if (isArray(v)) {
        if (v.length === 0) {
          newValue = undefined
        } else {
          newValue = v.map((item) => item instanceof Quantity ? item.toSI() : item)
        }
        k = `${k}:any`
      } else if (isPlainObject(v)) {
        newValue = cleanQuery(v, exclusive)
      } else {
        newValue = v
      }
    }
    newQuery[k] = newValue
  }
  return newQuery
}
