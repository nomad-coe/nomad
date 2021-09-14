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
import { useApi } from '../api'
import { setToArray, formatMeta, parseMeta } from '../../utils'
import searchQuantities from '../../searchQuantities'
import { Quantity } from '../../units'
import { useErrors } from '../errors'
import { combinePagination } from '../datatable/Datatable'

export const filters = new Set() // Contains the full names of all the available filters
export const filterGroups = [] // Mapping from filter full name -> group
export const filterAbbreviations = [] // Mapping of filter full name -> abbreviation
export const filterFullnames = [] // Mapping of filter abbreviation -> full name
export const filterData = {} // Stores data for each registered filter
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
 * This function is used to register a new filter within the FilterContext.
 * Filters are entities that can be searched throuh the filter panel and the
 * search bar, and can be encoded in the URL. Notice that a filter in this
 * context does not have to correspond to a quantity in the metainfo.
 *
 * Only registered filters may be searched for. The registration must happen
 * before any components use the filters. This is because:
 *  - The initial aggregation results must be fetched before any components
 *  using the filter values are rendered.
 *  - Several components need to know the list of available filters (e.g. the
 *  search bar and  the search panel). If filters are only registered during
 *  component initialization, it may already be too late to update other
 *  components.
 *
 * @param {string} name Name of the filter.
 * @param {string} group The group into which the filter belongs to. Groups
 * are used to e.g. in showing FilterSummaries about a group of filters.
 * @param {string|object} agg Custom setter/getter for the aggregation value. As a
 * shortcut you can provide an ES aggregation type as a string,
 * @param {object} value Custom setter/getter for the filter value.
 * @param {bool} multiple Whether this filter supports several values:
 * controls whether setting the value appends or overwrites.
 * @param {bool} exclusive Whether this filter is exclusive: only one value may be associated with an entry.
 */
function registerFilter(name, group, agg, value, multiple = true, exclusive = true) {
  filters.add(name)
  if (group) {
    filterGroups[group]
      ? filterGroups[group].add(name)
      : filterGroups[group] = new Set([name])
  }

  // Register mappings from full name to abbreviation and vice versa
  const abbreviation = name.split('.').pop()
  const oldName = filterAbbreviations[abbreviation]
  if (!oldName) {
    filterAbbreviations[name] = abbreviation
    filterFullnames[abbreviation] = name
  } else {
    delete filterFullnames[abbreviation]
    filterAbbreviations[name] = name
    filterAbbreviations[oldName] = oldName
  }

  const data = filterData[name] || {}
  if (agg) {
    let aggSet, aggGet
    if (isString(agg)) {
      aggSet = {[name]: agg}
      aggGet = (aggs) => (aggs[name][agg].data)
    } else {
      aggSet = agg.set
      aggGet = agg.get
    }
    data.aggSet = aggSet
    data.aggGet = aggGet
  }
  if (value) {
    data.valueSet = value.set
  }
  data.multiple = multiple
  data.exclusive = exclusive
  filterData[name] = data
}

// Filters that directly correspond to a metainfo value
registerFilter('results.material.structural_type', labelMaterial, 'terms')
registerFilter('results.material.functional_type', labelMaterial, 'terms')
registerFilter('results.material.compound_type', labelMaterial, 'terms')
registerFilter('results.material.material_name', labelMaterial)
registerFilter('results.material.chemical_formula_hill', labelElements)
registerFilter('results.material.chemical_formula_anonymous', labelElements)
registerFilter('results.material.n_elements', labelElements, 'min_max', undefined, false)
registerFilter('results.material.symmetry.bravais_lattice', labelSymmetry, 'terms')
registerFilter('results.material.symmetry.crystal_system', labelSymmetry, 'terms')
registerFilter('results.material.symmetry.structure_name', labelSymmetry, 'terms')
registerFilter('results.material.symmetry.strukturbericht_designation', labelSymmetry, 'terms')
registerFilter('results.material.symmetry.space_group_symbol', labelSymmetry)
registerFilter('results.material.symmetry.point_group', labelSymmetry)
registerFilter('results.material.symmetry.hall_symbol', labelSymmetry)
registerFilter('results.material.symmetry.prototype_aflow_id', labelSymmetry)
registerFilter('results.method.method_name', labelMethod, 'terms')
registerFilter('results.method.simulation.program_name', labelMethod, 'terms')
registerFilter('results.method.simulation.program_version', labelMethod)
registerFilter('results.method.simulation.dft.basis_set_type', labelDFT, 'terms')
registerFilter('results.method.simulation.dft.core_electron_treatment', labelDFT, 'terms')
registerFilter('results.method.simulation.dft.xc_functional_type', labelDFT, 'terms')
registerFilter('results.method.simulation.dft.relativity_method', labelDFT, 'terms')
registerFilter('results.method.simulation.gw.gw_type', labelGW, 'terms')
registerFilter('results.properties.electronic.band_structure_electronic.channel_info.band_gap_type', labelElectronic, 'terms')
registerFilter('results.properties.electronic.band_structure_electronic.channel_info.band_gap', labelElectronic, 'min_max', undefined, false)
registerFilter('external_db', labelAuthor, 'terms')
registerFilter('authors.name', labelAuthor)
registerFilter('upload_create_time', labelAuthor, 'min_max', undefined, false)
registerFilter('datasets.name', labelDataset)
registerFilter('datasets.doi', labelDataset)
registerFilter('entry_id', labelIDs)
registerFilter('upload_id', labelIDs)
registerFilter('results.material.material_id', labelIDs)
registerFilter('datasets.dataset_id', labelIDs)

// In exclusive element query the elements names are sorted and concatenated
// into a single string.
registerFilter(
  'results.material.elements',
  labelElements,
  'terms',
  {
    set: (newQuery, oldQuery, value) => {
      if (oldQuery.exclusive) {
        if (value.size !== 0) {
          newQuery['results.material.elements_exclusive'] = setToArray(value).sort().join(' ')
        }
      } else {
        newQuery['results.material.elements'] = value
      }
    }
  },
  true,
  false
)
// Electronic properties: subset of results.properties.available_properties
registerFilter(
  'electronic_properties',
  labelElectronic,
  {
    set: {'results.properties.available_properties': 'terms'},
    get: (aggs) => (aggs['results.properties.available_properties'].terms.data)
  },
  {
    set: (newQuery, oldQuery, value) => {
      const data = newQuery['results.properties.available_properties'] || new Set()
      value.forEach((item) => { data.add(item) })
      newQuery['results.properties.available_properties'] = data
    },
    get: (data) => (data.results.properties.available_properties)
  },
  true,
  false
)
// Vibrational properties: subset of results.properties.available_properties
registerFilter(
  'vibrational_properties',
  labelVibrational,
  {
    set: {'results.properties.available_properties': 'terms'},
    get: (aggs) => (aggs['results.properties.available_properties'].terms.data)
  },
  {
    set: (newQuery, oldQuery, value) => {
      const data = newQuery['results.properties.available_properties'] || new Set()
      value.forEach((item) => { data.add(item) })
      newQuery['results.properties.available_properties'] = data
    },
    get: (data) => (data.results.properties.available_properties)
  },
  true,
  false
)
// Visibility: controls the 'owner'-parameter in the API query, not part of the
// query itself.
registerFilter(
  'visibility',
  labelAccess,
  undefined,
  {set: () => {}},
  false
)
// Restricted: controls whether materials search is done in a restricted mode.
registerFilter(
  'restricted',
  undefined,
  undefined,
  {set: () => {}},
  false
)
// Exclusive: controls the way elements search is done.
registerFilter(
  'exclusive',
  undefined,
  undefined,
  {set: () => {}},
  false
)

// Material and entry queries target slightly different fields. Here we prebuild
// the mapping.
const materialNames = {} // Mapping of field name from entry -> material
const entryNames = {} // Mapping of field name from material -> entry
for (const name of Object.keys(searchQuantities)) {
  const prefix = 'results.material.'
  let materialName
  if (name.startsWith(prefix)) {
    materialName = name.substring(prefix.length)
  } else {
    materialName = `entries.${name}`
  }
  materialNames[name] = materialName
  entryNames[materialName] = name
}

export const searchContext = React.createContext()
export const SearchContext = React.memo(({
  resource,
  filtersLocked,
  children
}) => {
  const setQuery = useSetRecoilState(queryState)
  const setLocked = useSetRecoilState(lockedState)
  const {api} = useApi()
  const setInitialAggs = useSetRecoilState(initialAggsState)
  const [isMenuOpen, setIsMenuOpen] = useState(false)
  const [menuPath, setMenuPath] = useState('Filters')

  // Reset the query/locks when entering the search context for the first time
  const reset = useRecoilCallback(({reset}) => () => {
    for (let filter of filters) {
      reset(queryFamily(filter))
      reset(lockedFamily(filter))
    }
  }, [])

  useEffect(() => {
    reset()
  }, [reset])

  // Read the initial query from the URL
  const query = useMemo(() => {
    const location = window.location.href
    const split = location.split('?')
    let qs, query
    if (split.length === 1) {
      query = {}
    } else {
      qs = split.pop()
      query = qsToQuery(qs)
    }
    return query
  }, [])

  // Save the initial query and locked filters. Cannot be done inside useMemo
  // due to bad setState.
  useEffect(() => {
    setQuery(query)
    // Transform the locked values into a GUI-suitable format and store them
    if (filtersLocked) {
      const filtersLockedGUI = {}
      for (const [key, value] of Object.entries(filtersLocked)) {
        filtersLockedGUI[key] = toGUIFilter(key, value)
      }
      setLocked(filtersLockedGUI)
    }
  }, [setLocked, setQuery, query, filtersLocked])

  // Fetch initial aggregation data.
  useEffect(() => {
    const aggRequest = {}
    const aggNames = [...filters].filter(name => filterData[name].aggGet)
    for (const filter of aggNames) {
      toAPIAgg(aggRequest, filter, resource)
    }

    const search = {
      owner: 'visible',
      query: {},
      aggregations: aggRequest,
      pagination: {page_size: 0}
    }

    api.query(resource, search, false)
      .then(data => {
        const newData = {
          total: data.pagination.total,
          data: toGUIAgg(data.aggregations, aggNames, resource)
        }
        setInitialAggs(newData)
      })
  }, [api, setInitialAggs, resource])

  const values = useMemo(() => ({
    resource,
    isMenuOpen,
    setIsMenuOpen,
    menuPath,
    setMenuPath
  }), [resource, isMenuOpen, menuPath])

  return <searchContext.Provider value={values}>
    {children}
  </searchContext.Provider>
})
SearchContext.propTypes = {
  resource: PropTypes.string,
  filtersLocked: PropTypes.object,
  children: PropTypes.node
}

export function useSearchContext() {
  return useContext(searchContext)
}

/**
 * Each search filter is here mapped into a separate Recoil.js Atom. This
 * allows components to hook into individual search parameters (both for setting
 * and reading their value). This performs much better than having one large
 * Atom for the entire query, as this would cause all of the hooked components
 * to render even if they are not affected by some other search filter.
 * Re-renders became problematic with large and complex components (e.g. the
 * periodic table), for which the re-rendering takes significant time. Another
 * approach would have been to try and Memoize each sufficiently complex
 * component, but this quickly becomes a hard manual task.
 */
export const queryFamily = atomFamily({
  key: 'queryFamily',
  default: undefined
})
export const lockedFamily = atomFamily({
  key: 'lockedFamily',
  default: false
})

// Whether the search is initialized.
export const initializedState = atom({
  key: 'initialized',
  default: false
})

/**
 * Returns a function that can be called to reset all current filters.
 *
 * @returns Function for resetting all filters.
 */
export function useResetFilters() {
  const locked = useRecoilValue(lockedState)
  const reset = useRecoilCallback(({reset}) => () => {
    for (let filter of filters) {
      if (!locked[filter]) {
        reset(queryFamily(filter))
      }
    }
  }, [locked])
  return reset
}

/**
 * This hook will expose a function for reading if the given filter is locked.
 *
 * @param {string} name Name of the filter.
 * @returns Whether the filter is locked or not.
 */
export function useFilterLocked(name) {
  return useRecoilValue(lockedFamily(name))
}

/**
 * This hook will expose a function for reading the locked status of all
 * filters.
 *
 * @returns An object containing a mapping from filter name to a boolean
 * indicating whether it is locked or not.
 */
export function useFiltersLocked() {
  return useRecoilValue(lockedState)
}

/**
 * This hook will expose a function for reading if the given set of filters are
 * locked.
 *
 * @param {string} names Names of the filters.
 * @returns Array containing the filter values in a map and a setter function.
 */
let indexLocked = 0
export function useFiltersLockedState(names) {
  // We dynamically create a Recoil.js selector that is subscribed to the
  // filters specified in the input. This way only the specified filters will
  // cause a render. Recoil.js requires that each selector/atom has an unique
  // id. Because this hook can be called dynamically, we simply generate the ID
  // sequentially.
  const filterState = useMemo(() => {
    const id = `locked_selector${indexLocked}`
    indexLocked += 1
    return selector({
      key: id,
      get: ({get}) => {
        const query = {}
        for (let key of names) {
          const filter = get(lockedFamily(key))
          query[key] = filter
        }
        return query
      }
    })
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return useRecoilValue(filterState)
}

// Used to set the locked state of several filters at once
const lockedState = selector({
  key: 'lockedState',
  get: ({get}) => {
    const locks = {}
    for (let key of filters) {
      const filter = get(lockedFamily(key))
      locks[key] = filter
    }
    return locks
  },
  set: ({ get, set, reset }, data) => {
    if (data) {
      for (const [key, value] of Object.entries(data)) {
        set(queryFamily(key), value)
        set(lockedFamily(key), true)
      }
    }
  }
})

/**
 * This hook will expose a function for reading filter values. Use this hook if
 * you intend to only view the filter values and are not interested in setting
 * the filter.
 *
 * @param {string} name Name of the filter.
 * @returns currently set filter value.
 */
export function useFilterValue(name) {
  return useRecoilValue(queryFamily(name))
}

/**
 * This hook will expose a function for setting a filter value. Use this hook if
 * you intend to only set the filter value and are not interested in the query
 * results.
 *
 * @param {string} name Name of the quantity to set.
 * @returns function for setting the value for the given quantity
 */
export function useSetFilter(name) {
  return useSetRecoilState(queryFamily(name))
}

/**
 * This hook will expose a function for getting and setting filter values. Use
 * this hook if you intend to both read and write the filter value.
 *
 * @param {string} name Name of the filter.
 * @returns Array containing the filter value and setter function for it.
 */
export function useFilterState(name) {
  return useRecoilState(queryFamily(name))
}

/**
 * This hook will expose a function for setting the values of all filters.
 *
 * @returns An object containing a mapping from filter name to a boolean
 * indicating whether it is locked or not.
 */
export function useSetFilters() {
  return useSetRecoilState(filtersState)
}

// Used to get/set the locked state of all filters at once
const filtersState = selector({
  key: 'filtersState',
  get: ({get}) => {
    const query = {}
    for (let key of filters) {
      const filter = get(queryFamily(key))
      query[key] = filter
    }
    return query
  },
  set: ({set}, [key, value]) => {
    set(queryFamily(key), value)
  }
})

/**
 * This hook will expose a function for getting and setting filter values for
 * the specified list of filters. Use this hook if you intend to both read and
 * write the filter values.
 *
 * @param {string} names Names of the filters.
 * @returns Array containing the filter values in a map and a setter function.
 */
let indexFilters = 0
export function useFiltersState(names) {
  // We dynamically create a Recoil.js selector that is subscribed to the
  // filters specified in the input. This way only the specified filters will
  // cause a render. Recoil.js requires that each selector/atom has an unique
  // id. Because this hook can be called dynamically, we simply generate the ID
  // sequentially.
  const filterState = useMemo(() => {
    const id = `dynamic_selector${indexFilters}`
    indexFilters += 1
    return selector({
      key: id,
      get: ({get}) => {
        const query = {}
        for (let key of names) {
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
    for (let key of filters) {
      const filter = get(queryFamily(key))
      if (filter !== undefined) {
        query[key] = filter
      }
    }
    return query
  },
  set: ({ get, set, reset }, data) => {
    for (let filter of filters) {
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

  const updateQueryString = useCallback((query, locked) => {
    const queryString = queryToQs(query, locked)
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
  const query = qs.parse(queryString, {comma: true})
  const newQuery = {}
  for (let [key, value] of Object.entries(query)) {
    const split = key.split(':')
    key = split[0]
    let newKey = filterFullnames[key] || key
    const valueGUI = toGUIFilter(newKey, value)
    if (split.length !== 1) {
      const op = split[1]
      const oldValue = newQuery[newKey]
      if (!oldValue) {
        newQuery[newKey] = {[op]: valueGUI}
      } else {
        newQuery[newKey][op] = valueGUI
      }
    } else {
      newQuery[newKey] = valueGUI
    }
  }
  return newQuery
}

/**
 * Converts a query into an object that uses NOMAD API query string keys and values.
 * @param {object} query A query object representing the currently active
 * filters.
 * @returns {object} An object that can be serialized into a query string
 */
export function queryToQsData(query, locked) {
  locked = locked || {}
  const queryStringQuery = {}
  for (const [key, value] of Object.entries(query)) {
    if (locked[key]) {
      continue
    }
    const {formatter} = formatMeta(key, false)
    let newValue
    const newKey = filterAbbreviations[key] || key
    if (isPlainObject(value)) {
      if (!isNil(value.gte)) {
        queryStringQuery[`${newKey}:gte`] = formatter(value.gte)
      }
      if (!isNil(value.lte)) {
        queryStringQuery[`${newKey}:lte`] = formatter(value.lte)
      }
    } else {
      if (isArray(value)) {
        newValue = value.map(formatter)
      } else if (value instanceof Set) {
        newValue = [...value].map(formatter)
      } else {
        newValue = formatter(value)
      }
      queryStringQuery[newKey] = newValue
    }
  }
  return queryStringQuery
}

/**
 * Converts a query into a valid query string.
 * @param {object} query A query object representing the currently active
 * filters.
 * @returns URL querystring, not encoded if possible to improve readability.
 */
function queryToQs(query, locked) {
  const queryData = queryToQsData(query, locked)
  return qs.stringify(queryData, {indices: false, encode: false})
}

export const initialAggsState = atom({
  key: 'initialAggs',
  default: undefined
})

/**
 * Hook for returning an initial aggregation value for a filter.
 *
 * @returns {array} Array containing the aggregation data.
 */
export function useInitialAgg(name) {
  const aggs = useRecoilValue(initialAggsState)
  return aggs?.[name]
}

/**
 * Hook for retrieving the most up-to-date aggregation results for a specific
 * filter, taking into account the current search context.
 *
 * @param {string} name The filter name
 * @param {bool} update Whether the hook needs to react to changes in the
 * current query context. E.g. if the component showing the data is not visible,
 * this can be set to false.
 *
 * @returns {array} The data-array returned by the API.
 */
export function useAgg(name, update = true, delay = 500) {
  const {api} = useApi()
  const { resource } = useSearchContext()
  const [results, setResults] = useState(undefined)
  const initialAggs = useRecoilValue(initialAggsState)
  const restrict = filterData[name].exclusive
  const query = useQuery()
  const firstLoad = useRef(true)

  // Pretty much all of the required pre-processing etc. should be done in this
  // function, as it is the final one that gets called after the debounce
  // interval.
  const apiCall = useCallback((query) => {
    // If the restrict option is enabled, the filters targeting the specified
    // quantity will be removed. This way all possible options pre-selection can
    // be returned.
    let queryCleaned = {...query}
    if (restrict && query && name in query) {
      delete queryCleaned[name]
    }
    queryCleaned = toAPIQuery(queryCleaned, resource, query.restricted)
    const aggRequest = {}
    toAPIAgg(aggRequest, name, resource)
    const search = {
      owner: query.visibility || 'visible',
      query: queryCleaned,
      aggregations: aggRequest,
      pagination: {page_size: 0},
      required: { include: [] }
    }

    api.query(resource, search, false)
      .then(data => {
        const newData = {
          total: data.pagination.total,
          data: toGUIAgg(data.aggregations, [name], resource)[name]
        }
        firstLoad.current = false
        setResults(newData)
      })
  }, [api, name, restrict, resource])

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
      if (isEmpty(query) && !isNil(initialAggs?.[name])) {
        setResults(initialAggs[name])
      // Make an immediate request for the aggregation values if query has been
      // specified.
      } else {
        apiCall(query)
      }
    } else {
      debounced(query)
    }
  }, [apiCall, name, debounced, query, update, initialAggs])

  return results
}

/**
 * Hook that adds results and pagination state to your component.
 *
 * @param {number} delay The debounce delay in milliseconds.
 *
 * @returns {object} An object with data (the search results as array), pagination
 *  (the current api pagination object), and setPagination (to update pagination)
 */
export function useScrollResults(initialPagination, delay = 500) {
  const {api} = useApi()
  const {raiseErrors} = useErrors()
  const {resource} = useSearchContext()
  const query = useQuery()
  const locked = useRecoilValue(lockedState)
  const updateQueryString = useUpdateQueryString()

  const [results, setResults] = useState([])
  const paginationResponse = useRef({})
  const [pagination, setPagination] = useState(initialPagination || {page_size: 10})

  // This callback will call the API with a debounce by the given delay.
  const callApi = useCallback(debounce((query, locked, pagination) => {
    const isExtend = pagination.page_after_value

    const request = {
      owner: query.visibility || 'visible',
      query: toAPIQuery(query, resource, query.restricted),
      pagination: pagination
    }

    api.query(resource, request)
      .then(response => {
        paginationResponse.current = response.pagination
        if (isExtend) {
          setResults(results => [...results, ...response.data])
        } else {
          setResults(response.data)
        }

        // We only update the query string after the API call is finished. Updating
        // the query string causes quite an intensive render (not sure why), so it
        // is better to debounce this value as well to keep the user interaction
        // smoother.
        updateQueryString(query, locked)
      })
      .catch(raiseErrors)
  }, delay), [resource, api, raiseErrors, updateQueryString, paginationResponse])

  // Whenever the query or pagination changes, we make an api call.
  // The results are fetched as a side effect in order to not block the
  // rendering. This causes two renders: first one without the data, the second
  // one with the data.
  useEffect(() => {
    if (!query) {
      return
    }

    if (pagination.page_after_value && pagination.page_after_value === paginationResponse.current.page_after_value) {
      pagination.page_after_value = undefined
      pagination.next_page_after_value = undefined
    }

    callApi(query, locked, pagination)
  }, [callApi, query, locked, pagination])

  return {
    data: results,
    pagination: combinePagination(pagination, paginationResponse.current),
    setPagination: setPagination}
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
export function toAPIQuery(query, resource, restricted) {
  // Perform custom transformations
  let queryCustomized = {}
  for (let [k, v] of Object.entries(query)) {
    const setter = filterData[k]?.valueSet
    if (setter) {
      setter(queryCustomized, query, v)
    } else {
      queryCustomized[k] = v
    }
  }

  let queryNormalized = {}
  for (const [k, v] of Object.entries(queryCustomized)) {
    // Transform sets into lists and Quantities into SI values and modify keys
    // according to target resource (entries/materials).
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

    // The postfixes are added here. By default query items with array values
    // get the 'any'-postfix.
    let postfix
    if (isArray(newValue)) {
      const fieldPostfixMap = {
        'results.properties.available_properties': 'all',
        'results.material.elements': 'all'
      }
      postfix = fieldPostfixMap[k] || 'any'
    }

    // For material query the keys are remapped.
    let newKey = resource === 'materials' ? materialNames[k] : k
    newKey = postfix ? `${newKey}:${postfix}` : newKey
    queryNormalized[newKey] = newValue
  }

  if (resource === 'materials') {
    // In restricted search we simply move all method/properties filters
    // inside a single entries-subsection.
    if (restricted) {
      const entrySearch = {}
      for (const [k, v] of Object.entries(queryNormalized)) {
        if (k.startsWith('entries.')) {
          const name = k.split('entries.').pop()
          entrySearch[name] = v
          delete queryNormalized[k]
        }
      }
      if (!isEmpty(entrySearch)) {
        queryNormalized.entries = entrySearch
      }
    // In unrestricted search we have to split each filter and each filter value
    // into it's own separate entries query. These queries are then joined with
    // 'and'.
    } else {
      const entrySearch = []
      for (const [k, v] of Object.entries(queryNormalized)) {
        if (k.startsWith('entries.')) {
          const newKey = k.split(':')[0]
          if (isArray(v)) {
            for (const item of v) {
              entrySearch.push({[newKey]: item})
            }
          } else {
            entrySearch.push({[newKey]: v})
          }
          delete queryNormalized[k]
        }
      }
      if (entrySearch.length > 0) {
        queryNormalized.and = entrySearch
      }
    }
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
 * Cleans a filter value into a form that is supported by the GUI. This includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 *
 * @returns {any} The filter value in a format that is suitable for the GUI.
 */
export function toGUIFilter(name, value, units = undefined) {
  let multiple = filterData[name].multiple
  let newValue
  const {parser} = parseMeta(name)
  if (isArray(value)) {
    newValue = new Set(value.map((v) => parser(v, units)))
  } else if (isPlainObject(value)) {
    newValue = {}
    if (!isNil(value.gte)) {
      newValue.gte = parser(value.gte, units)
    }
    if (!isNil(value.lte)) {
      newValue.lte = parser(value.lte, units)
    }
  } else {
    newValue = parser(value, units)
    if (multiple) {
      newValue = new Set([newValue])
    }
  }
  return newValue
}

/**
 * Used to transform a GUI aggregation query into a form that is usable by the
 * API.
 *
 * @param {object} aggs The aggregation data in which the modifications are
 * made.
 * @param {string} filter The filter name
 * @param {string} resource The resource we are looking at: entries or materials.
 */
function toAPIAgg(aggs, filter, resource) {
  const aggSet = filterData[filter].aggSet
  if (aggSet) {
    for (const [key, type] of Object.entries(aggSet)) {
      const name = resource === 'materials' ? materialNames[key.split(':')[0]] : key
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
 * @param {array} filters The filters to take into account.
 * @param {string} resource The resource we are looking at: entries or materials.
 *
 * @returns {object} Aggregation data that is usable by the GUI.
 */
function toGUIAgg(aggs, filters, resource) {
  if (isEmpty(aggs)) {
    return aggs
  }
  // Modify keys according to target resource (entries/materials).
  let aggsNormalized
  if (resource === 'materials') {
    aggsNormalized = {}
    for (const key of Object.keys(aggs)) {
      const name = resource === 'materials' ? entryNames[key] : key
      aggs[key].quantity = name
      aggsNormalized[name] = aggs[key]
    }
  } else {
    aggsNormalized = aggs
  }

  // Perform custom transformations
  const aggsCustomized = {}
  for (const name of filters) {
    const aggGet = filterData[name].aggGet
    if (aggGet) {
      let agg
      agg = aggGet(aggsNormalized)
      aggsCustomized[name] = agg
    }
  }
  return aggsCustomized
}
