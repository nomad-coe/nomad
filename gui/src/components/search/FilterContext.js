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
import React, { useCallback, useEffect, useState, useRef, useMemo, useContext } from 'react'
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
  isNil,
  isString
} from 'lodash'
import qs from 'qs'
import PropTypes from 'prop-types'
import { useHistory } from 'react-router-dom'
import { useApi } from '../apiV1'
import { setToArray, formatMeta, parseMeta } from '../../utils'
import searchQuantities from '../../searchQuantities'
import { Quantity } from '../../units'

let index = 0
export const quantities = new Set()
export const quantityGroups = new Map()
export const quantityAggregations = {}
export const quantityAggGets = {}
export const quantityAggSets = {}
export const quantityAbbreviations = new Map()
export const quantityFullnames = new Map()
export const quantityMaterialNames = {}
export const quantityEntryNames = {}
export const quantityData = {}
export const labelMaterial = 'Material'
export const labelElements = 'Elements / Formula'
export const labelSymmetry = 'Symmetry'
export const labelMethod = 'Method'
export const labelSimulation = 'Simulation'
export const labelDFT = 'DFT'
export const labelGW = 'GW'
export const labelProperties = 'Properties'
export const labelElectronic = 'Electronic'
export const labelVibrational = 'Vibrational'
export const labelAuthor = 'Author / Origin'
export const labelAccess = 'Access'
export const labelDataset = 'Dataset'
export const labelIDs = 'IDs'

/**
 * This function is used to register a "quantity" within the FilterContext.
 * Quantities are entities that can be searched throuh the filter panel and the
 * search bar, and can be encoded in the URL. Notice that a quantity in this
 * context does not have to correspond to a quantity in the metainfo.
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
 * @param {string|object} agg Custom setter/getter for the aggregation value. As a
 * shortcut you can provide an ES aggregation type as a string,
 * @param {object} value Custom setter/getter for the quantity value.
 * @param {bool} multiple Whether this quantity supports several values:
 * controls whether setting the value appends or overwrites.
 */
function registerQuantity(name, group, agg, value, multiple = true) {
  // Store the available quantities, their grouping, and the initial aggregation
  // types.
  quantities.add(name)
  quantityGroups.has(group) ? quantityGroups.get(group).add(name) : quantityGroups.set(group, new Set([name]))

  // Register mappings from full name to abbreviation and vice versa
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

  // Store the string/function that is used to getch the aggregation data for
  // this quantity.
  if (agg) {
    let aggSet, aggGet
    if (isString(agg)) {
      aggSet = {[name]: agg}
      aggGet = (aggs) => (aggs[name][agg].data)
    } else {
      aggSet = agg.set
      aggGet = agg.get
    }
    quantityAggSets[name] = aggSet
    quantityAggGets[name] = aggGet
  }

  // Store the default query postfixes (any/all) for each quantity
  const data = quantityData[name] || {}
  if (value) {
    data.valueGet = value.get
    data.valueSet = value.set
  }
  data.multiple = multiple
  quantityData[name] = data
}

// Quantities that directly correspond to a metainfo value
registerQuantity('results.material.structural_type', labelMaterial, 'terms')
registerQuantity('results.material.functional_type', labelMaterial, 'terms')
registerQuantity('results.material.compound_type', labelMaterial, 'terms')
registerQuantity('results.material.material_name', labelMaterial)
registerQuantity('results.material.chemical_formula_hill', labelElements)
registerQuantity('results.material.chemical_formula_anonymous', labelElements)
registerQuantity('results.material.n_elements', labelElements, 'min_max', undefined, false)
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
registerQuantity('results.method.simulation.dft.xc_functional_type', labelDFT, 'terms')
registerQuantity('results.method.simulation.dft.relativity_method', labelDFT, 'terms')
registerQuantity('results.method.simulation.gw.gw_type', labelGW, 'terms')
registerQuantity('results.properties.electronic.band_structure_electronic.channel_info.band_gap_type', labelElectronic, 'terms')
registerQuantity('results.properties.electronic.band_structure_electronic.channel_info.band_gap', labelElectronic, 'min_max', undefined, false)
registerQuantity('external_db', labelAuthor, 'terms')
registerQuantity('authors.name', labelAuthor)
registerQuantity('upload_time', labelAuthor, 'min_max', undefined, false)
registerQuantity('datasets.name', labelDataset)
registerQuantity('datasets.doi', labelDataset)
registerQuantity('entry_id', labelIDs)
registerQuantity('upload_id', labelIDs)
registerQuantity('results.material.material_id', labelIDs)
registerQuantity('datasets.dataset_id', labelIDs)

// In regular element query we use the 'all'-postfix.
registerQuantity(
  'results.material.elements',
  labelElements,
  'terms',
  { set: (query, value) => (query['results.material.elements:all'] = value) }
)
// In exclusive element query the elements names are sorted and concatenated
// into a single string.
registerQuantity(
  'results.material.elements_exclusive',
  labelElements,
  undefined,
  {
    set: (query, value) => {
      if (value.size !== 0) {
        query['results.material.elements_exclusive'] = setToArray(value).sort().join(' ')
      }
    }
  }
)
// Electronic properties: subset of results.properties.available_properties
registerQuantity(
  'electronic_properties',
  labelElectronic,
  {
    set: {'results.properties.available_properties': 'terms'},
    get: (aggs) => (aggs['results.properties.available_properties'].terms.data)
  },
  {
    set: (query, value) => {
      const data = query['results.properties.available_properties'] || new Set()
      value.forEach((item) => { data.add(item) })
      query['results.properties.available_properties'] = data
    },
    get: (data) => (data.results.properties.available_properties)
  }
)
// Vibrational properties: subset of results.properties.available_properties
registerQuantity(
  'vibrational_properties',
  labelVibrational,
  {
    set: {'results.properties.available_properties': 'terms'},
    get: (aggs) => (aggs['results.properties.available_properties'].terms.data)
  },
  {
    set: (query, value) => {
      const data = query['results.properties.available_properties'] || new Set()
      value.forEach((item) => { data.add(item) })
      query['results.properties.available_properties'] = data
    },
    get: (data) => (data.results.properties.available_properties)
  }
)
// Visibility: controls the 'owner'-parameter in the API query, not part of the
// query itself.
registerQuantity(
  'visibility',
  labelAccess,
  undefined,
  {
    set: () => {},
    get: () => {}
  },
  false
)

// Material and entry queries target slightly different fields. Here we prebuild
// the mapping.
for (const name of Object.keys(searchQuantities)) {
  const prefix = 'results.material.'
  let materialName
  if (name.startsWith(prefix)) {
    materialName = name.substring(prefix.length)
  } else {
    materialName = `entries.${name}`
  }
  quantityMaterialNames[name] = materialName
  quantityEntryNames[materialName] = name
}

export const searchContext = React.createContext()
export const SearchContext = React.memo(({
  children
}) => {
  const setQuery = useSetRecoilState(queryState)
  const {api} = useApi()
  const setInitialAggs = useSetRecoilState(initialAggsState)

  // Reset the query when entering the search context for the first time
  useEffect(() => {
    setQuery(undefined)
  }, [setQuery])

  // Read the target resource and initial query from the URL
  const [resource, query] = useMemo(() => {
    const location = window.location.href
    const split = location.split('?')
    let qs, path, query
    if (split.length === 1) {
      path = split.pop()
      query = {}
    } else {
      [path, qs] = split
      query = qsToQuery(qs)
    }
    return [path.split('/').pop(), query]
  }, [])

  // Save the initial query. Cannot be done inside useMemo due to bad setState.
  useEffect(() => {
    setQuery(query)
  }, [setQuery, query])

  // Fetch initial aggregation data. We include all data here.
  useEffect(() => {
    const aggRequest = {}
    const aggNames = Object.keys(quantityAggGets)
    for (const quantity of aggNames) {
      toAPIAgg(aggRequest, quantity, resource)
    }

    const search = {
      owner: 'visible',
      query: {},
      aggregations: aggRequest,
      pagination: {page_size: 0}
    }

    api.query(resource, search, false)
      .then(data => {
        data = toGUIAgg(data.aggregations, aggNames, resource)
        setInitialAggs(data)
      })
  }, [api, setInitialAggs, resource])

  const values = useMemo(() => ({
    resource
  }), [resource])

  return <searchContext.Provider value={values}>
    {children}
  </searchContext.Provider>
})
SearchContext.propTypes = {
  children: PropTypes.node
}

export function useSearchContext() {
  return useContext(searchContext)
}

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

// Current menu path
export const menuPath = atom({
  key: 'menuPath',
  default: 'Filters'
})
export function useMenuPathState() {
  return useRecoilState(menuPath)
}
export function useMenuPath() {
  return useRecoilValue(menuPath)
}
export function useSetMenuPath() {
  return useSetRecoilState(menuPath)
}

// Whether the search is initialized.
export const initializedState = atom({
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
    if (!get(initializedState)) {
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
    for (let filter of quantities) {
      reset(queryFamily(filter))
    }
    if (data) {
      for (const [key, value] of Object.entries(data)) {
        set(queryFamily(key), value)
      }
      set(initializedState, true)
    } else {
      set(initializedState, false)
    }
  }
})

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
        if (type !== 'number' && type !== 'timestamp' && key !== 'visibility') {
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

export const initialAggsState = atom({
  key: 'initialAggs',
  default: undefined
})

/**
 * Hook for returning the current search object.
 *
 * @returns {object} Object containing the search object.
 */
export function useInitialAgg(quantity) {
  const aggs = useRecoilValue(initialAggsState)
  return aggs?.[quantity]
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
export function useAgg(quantity, restrict = false, update = true, delay = 500) {
  const {api} = useApi()
  const { resource } = useSearchContext()
  const [results, setResults] = useState(undefined)
  const initialAggs = useRecoilValue(initialAggsState)
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
      delete queryCleaned[quantity]
    }
    queryCleaned = toAPIQuery(queryCleaned, resource)
    const aggRequest = {}
    toAPIAgg(aggRequest, quantity, resource)
    const search = {
      owner: query.visibility || 'visible',
      query: queryCleaned,
      aggregations: aggRequest,
      pagination: {page_size: 0},
      required: { include: [] }
    }

    api.query(resource, search, false)
      .then(data => {
        data = toGUIAgg(data.aggregations, [quantity], resource)
        firstLoad.current = false
        setResults(data[quantity])
      })
  }, [api, quantity, restrict, resource])

  // This is a debounced version of apiCall.
  const debounced = useCallback(debounce(apiCall, delay), [])

  // The API call is made immediately on first render. On subsequent renders it
  // will be debounced.
  useEffect(() => {
    if (!update || query === undefined) {
      return
    }
    if (firstLoad.current) {
      // Fetch the initial aggregation values if no query
      // is specified.
      if (isEmpty(query)) {
        setResults(initialAggs[quantity])
      // Make an immediate request for the aggregation values if query has been
      // specified.
      } else {
        apiCall(query, exclusive)
      }
    } else {
      debounced(query, exclusive)
    }
  }, [apiCall, quantity, debounced, query, exclusive, update, initialAggs])

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
export function useScrollResults(pageSize, orderBy, order, exclusive, delay = 500) {
  const {api} = useApi()
  const {resource} = useSearchContext()
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
    const cleanedQuery = toAPIQuery(query, resource)
    const search = {
      owner: query.visibility || 'visible',
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
    api.query(resource, search)
      .then(data => {
        pageAfterValue.current = data.pagination.next_page_after_value
        total.current = data.pagination.total
        setResults(data)
        loading.current = false
      })

    // We only update the query string after the API call is finished. Updating
    // the query string causes quite an intensive render (not sure why), so it
    // is better to debounce this value as well to keep the user interaction
    // smoother.
    updateQueryString(query)
  }, [api, updateQueryString, resource])

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
    api.query(resource, searchRef.current)
      .then(data => {
        pageAfterValue.current = data.pagination.next_page_after_value
        total.current = data.pagination.total
        setResults(old => {
          data.data = old.data.concat(data.data)
          return data
        })
        loading.current = false
      })
  }, [api, resource])

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
      debounced(query, pageSize, orderBy, order, exclusive)
    }
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
 * Converts the contents a query into a format that is suitable for the API.
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
export function toAPIQuery(query, resource) {
  // Perform custom transformations
  let queryCustomized = {}
  for (let [k, v] of Object.entries(query)) {
    const setter = quantityData[k]?.valueSet
    if (setter) {
      setter(queryCustomized, v)
    } else {
      queryCustomized[k] = v
    }
  }

  // Transform sets into lists and Quantities into SI values and modify keys
  // according to target resource (entries/materials).
  let queryNormalized = {}
  for (let [k, v] of Object.entries(queryCustomized)) {
    let newValue
    if (isPlainObject(v)) {
      newValue = {}
      if (!isNil(v.lte)) {
        newValue.lte = toAPIQueryValue(v.lte)
      }
      if (!isNil(v.gte)) {
        newValue.gte = toAPIQueryValue(v.gte)
      }
    } else {
      newValue = toAPIQueryValue(v)
    }
    k = resource === 'materials' ? quantityMaterialNames[k.split(':')[0]] : k
    queryNormalized[k] = newValue
  }

  return queryNormalized
}

/**
 * Cleans a filter value into a form that is supported by the API. This includes:
 * - Sets are transformed into Arrays
 * - Quantities are converted to SI values.
 *
 * @returns {any} The filter value in a format that is suitable for the API.
 */
function toAPIQueryValue(value) {
  let newValue
  if (value instanceof Set) {
    newValue = setToArray(value)
    if (newValue.length === 0) {
      newValue = undefined
    } else {
      newValue = newValue.map((item) => item instanceof Quantity ? item.toSI() : item)
    }
  } else if (value instanceof Quantity) {
    newValue = value.toSI()
  } else if (isArray(value)) {
    if (value.length === 0) {
      newValue = undefined
    } else {
      newValue = value.map((item) => item instanceof Quantity ? item.toSI() : item)
    }
  } else {
    newValue = value
  }
  return newValue
}

/**
 * Used to transform a GUI aggregation query into a form that is usable by the
 * API.
 *
 * @param {object} aggs The aggregation data in which the modifications are
 * made.
 * @param {string} quantity The quantity name
 * @param {string} resource The resource we are looking at: entries or materials.
 */
function toAPIAgg(aggs, quantity, resource) {
  const aggSet = quantityAggSets[quantity]
  if (aggSet) {
    for (const [key, type] of Object.entries(aggSet)) {
      const name = resource === 'materials' ? quantityMaterialNames[key.split(':')[0]] : key
      const agg = aggs[name] || {}
      agg[type] = {
        quantity: name,
        size: 500
      }
      aggs[name] = agg
    }
  }
}

/**
 * Used to transform an API aggregation query into a form that is usable by the
 * GUI.
 *
 * @param {object} aggs The aggregation data as returned by the API.
 * @param {string} quantity The quantity name
 * @param {string} resource The resource we are looking at: entries or materials.
 *
 * @returns {object} Aggregation data for the given quantity.
 */
function toGUIAgg(aggs, quantities, resource) {
  if (isEmpty(aggs)) {
    return aggs
  }
  // Modify keys according to target resource (entries/materials).
  let aggsNormalized
  if (resource === 'materials') {
    aggsNormalized = {}
    for (const key of Object.keys(aggs)) {
      const name = resource === 'materials' ? quantityEntryNames[key] : key
      aggs[key].quantity = name
      aggsNormalized[name] = aggs[key]
    }
  } else {
    aggsNormalized = aggs
  }

  // Perform custom transformations
  const aggsCustomized = {}
  for (const name of quantities) {
    const aggGet = quantityAggGets[name]
    if (aggGet) {
      let agg
      agg = aggGet(aggsNormalized)
      aggsCustomized[name] = agg
    }
  }
  return aggsCustomized
}
