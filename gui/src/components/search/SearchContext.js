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
import InputList from './input/InputList'
import InputSlider from './input/InputSlider'
import InputDateRange from './input/InputDateRange'
import InputPeriodicTable from './input/InputPeriodicTable'

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
function registerFilter(name, group, statConfig, agg, value, multiple = true, exclusive = true, options) {
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
  data.statConfig = statConfig
  data.options = options
  filterData[name] = data
}

const listStatConfig = {
  component: InputList,
  layout: {
    widthOptions: ['small', 'medium', 'large'],
    widthDefault: 'small',
    ratioDefault: 3 / 4,
    ratioOptions: [1, 3 / 4]
  }
}
const ptStatConfig = {
  component: InputPeriodicTable,
  layout: {
    widthOptions: ['small', 'medium', 'large'],
    widthDefault: 'large',
    ratioOptions: [3 / 2],
    ratioDefault: 3 / 2
  }
}
export const widthMapping = {
  large: {
    'sm': 12,
    'md': 9,
    'lg': 8,
    'xl': 6
  },
  medium: {
    'sm': 6,
    'md': 6,
    'lg': 4,
    'xl': 3
  },
  small: {
    'sm': 6,
    'md': 3,
    'lg': 4,
    'xl': 3
  }
}

// Filters that directly correspond to a metainfo value
registerFilter('results.material.structural_type', labelMaterial, listStatConfig, 'terms')
registerFilter('results.material.functional_type', labelMaterial, listStatConfig, 'terms')
registerFilter('results.material.compound_type', labelMaterial, listStatConfig, 'terms')
registerFilter('results.material.material_name', labelMaterial, listStatConfig)
registerFilter('results.material.chemical_formula_hill', labelElements, listStatConfig)
registerFilter('results.material.chemical_formula_anonymous', labelElements, listStatConfig)
registerFilter('results.material.n_elements', labelElements, InputSlider, 'min_max', undefined, false)
registerFilter('results.material.symmetry.bravais_lattice', labelSymmetry, listStatConfig, 'terms')
registerFilter('results.material.symmetry.crystal_system', labelSymmetry, listStatConfig, 'terms')
registerFilter('results.material.symmetry.structure_name', labelSymmetry, listStatConfig, 'terms')
registerFilter('results.material.symmetry.strukturbericht_designation', labelSymmetry, listStatConfig, 'terms')
registerFilter('results.material.symmetry.space_group_symbol', labelSymmetry, listStatConfig)
registerFilter('results.material.symmetry.point_group', labelSymmetry, listStatConfig)
registerFilter('results.material.symmetry.hall_symbol', labelSymmetry, listStatConfig)
registerFilter('results.material.symmetry.prototype_aflow_id', labelSymmetry, listStatConfig)
registerFilter('results.method.method_name', labelMethod, listStatConfig, 'terms')
registerFilter('results.method.simulation.program_name', labelMethod, listStatConfig, 'terms')
registerFilter('results.method.simulation.program_version', labelMethod, listStatConfig)
registerFilter('results.method.simulation.dft.basis_set_type', labelDFT, listStatConfig, 'terms')
registerFilter('results.method.simulation.dft.core_electron_treatment', labelDFT, listStatConfig, 'terms')
registerFilter('results.method.simulation.dft.xc_functional_type', labelDFT, listStatConfig, 'terms')
registerFilter('results.method.simulation.dft.relativity_method', labelDFT, listStatConfig, 'terms')
registerFilter('results.method.simulation.gw.gw_type', labelGW, listStatConfig, 'terms')
registerFilter('results.properties.electronic.band_structure_electronic.channel_info.band_gap_type', labelElectronic, listStatConfig, 'terms')
registerFilter('results.properties.electronic.band_structure_electronic.channel_info.band_gap', labelElectronic, InputSlider, 'min_max', undefined, false)
registerFilter('external_db', labelAuthor, listStatConfig, 'terms')
registerFilter('authors.name', labelAuthor, listStatConfig)
registerFilter('upload_create_time', labelAuthor, InputDateRange, 'min_max', undefined, false)
registerFilter('datasets.name', labelDataset, listStatConfig)
registerFilter('datasets.doi', labelDataset, listStatConfig)
registerFilter('entry_id', labelIDs, listStatConfig)
registerFilter('upload_id', labelIDs, listStatConfig)
registerFilter('results.material.material_id', labelIDs, listStatConfig)
registerFilter('datasets.dataset_id', labelIDs, listStatConfig)

// In exclusive element query the elements names are sorted and concatenated
// into a single string.
registerFilter(
  'results.material.elements',
  labelElements,
  ptStatConfig,
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
const electronicOptions = {
  band_structure_electronic: {label: 'band structure'},
  dos_electronic: {label: 'density of states'}
}
const electronicProps = new Set(Object.keys(electronicOptions))
registerFilter(
  'electronic_properties',
  labelElectronic,
  listStatConfig,
  {
    set: {'results.properties.available_properties': 'terms'},
    get: (aggs) => (aggs['results.properties.available_properties'].terms.data
      .filter((value) => electronicProps.has(value.value)))
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
  false,
  electronicOptions
)
// Vibrational properties: subset of results.properties.available_properties
export const vibrationalOptions = {
  dos_phonon: {label: 'phonon density of states'},
  band_structure_phonon: {label: 'phonon band structure'},
  energy_free_helmholtz: {label: 'Helmholtz free energy'},
  heat_capacity_constant_volume: {label: 'heat capacity constant volume'}
}
const vibrationalProps = new Set(Object.keys(vibrationalOptions))
registerFilter(
  'vibrational_properties',
  labelVibrational,
  listStatConfig,
  {
    set: {'results.properties.available_properties': 'terms'},
    get: (aggs) => (aggs['results.properties.available_properties'].terms.data
      .filter((value) => vibrationalProps.has(value.value)))
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
  false,
  vibrationalOptions
)
// Visibility: controls the 'owner'-parameter in the API query, not part of the
// query itself.
registerFilter(
  'visibility',
  labelAccess,
  undefined,
  undefined,
  {set: () => {}},
  false
)
// Restricted: controls whether materials search is done in a restricted mode.
registerFilter(
  'restricted',
  undefined,
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

/**
 * React context that provides access to the search state implemented with
 * Recoil.js. The purpose of this Context is to hide the Recoil.js
 * implementation details and provide a clean access to individual states in
 * order to prevent unnecessary re-renders (if we were to provide state values
 * and setters in a vanilla React Context, updating a single filter would cause
 * all components using the context to update.)
 *
 * At the moment we use a global set of Recoil.js atoms. These are reused throughout the application.
 * If in the future we need to use several searchcontexts simultaneously, each
 * SearchContext could instantiate it's own set of atoms dynamically.
 */
export const searchContext = React.createContext()
export const SearchContext = React.memo(({
  resource,
  filtersLocked,
  children
}) => {
  const setQuery = useSetRecoilState(queryState)
  const setLocked = useSetRecoilState(lockedState)
  const setStatistics = useSetRecoilState(statisticsState)
  const {api} = useApi()
  const setInitialAggs = useSetRecoilState(initialAggsState)

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
  const [query, statistics] = useMemo(() => {
    const location = window.location.href
    const split = location.split('?')
    let qs, query, statistics
    if (split.length === 1) {
      query = {}
    } else {
      qs = split.pop();
      [query, statistics] = qsToSearch(qs)
      console.log(statistics)
    }
    return [query, statistics]
  }, [])

  // Save the initial query, locked filters and statistics. Cannot be done
  // inside useMemo due to bad setState.
  useEffect(() => {
    setQuery(query)
    setStatistics(statistics)
    // Transform the locked values into a GUI-suitable format and store them
    if (filtersLocked) {
      const filtersLockedGUI = {}
      for (const [key, value] of Object.entries(filtersLocked)) {
        filtersLockedGUI[key] = toGUIFilter(key, value)
      }
      setLocked(filtersLockedGUI)
    }
  }, [setLocked, setQuery, setStatistics, query, statistics, filtersLocked])

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
    useIsMenuOpen: () => useRecoilValue(isMenuOpenState),
    useSetIsMenuOpen: () => useSetRecoilState(isMenuOpenState),
    useIsStatisticsEnabled: () => useRecoilValue(isStatisticsEnabledState),
    useSetIsStatisticsEnabled: () => useSetRecoilState(isStatisticsEnabledState)
  }), [resource])

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
export const initializedState = atom({
  key: 'initialized',
  default: false
})
export const isStatisticsEnabledState = atom({
  key: 'statisticsEnabled',
  default: true
})
export const isMenuOpenState = atom({
  key: 'isMenuOpen',
  default: false
})

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

export const statisticFamily = atomFamily({
  key: 'statisticFamily',
  default: undefined
})

// Used to get/set the the statistics configuration of all filters
const statisticsState = selector({
  key: 'statisticsState',
  set: ({set}, stats) => {
    if (stats) {
      for (let [key, value] of Object.entries(stats)) {
        set(statisticFamily(key), value)
      }
    }
  },
  get: ({get}) => {
    const stats = {}
    for (let filter of filters) {
      const stat = get(statisticFamily(filter))
      if (stat) stats[filter] = stat
    }
    return stats
  }
})

/**
 * This hook will expose a function for reading if filter statistics are shown
 * on the search page. Use this hook if you intend to only view the value and
 * are not interested in setting the value.
 *
 * @param {string} name Name of the filter.
 * @returns Whether the filter statistics are shown.
 */
export function useStatisticValue(name) {
  return useRecoilValue(statisticFamily(name))
}

/**
 * This hook will expose a function for setting if filter statistics are shown
 * on the search page. Use this hook if you intend to only set the value and are
 * not interested in reading it.
 *
 * @param {string} name Name of the quantity to set.
 * @returns function for setting the value
 */
export function useSetStatistic(name) {
  return useSetRecoilState(statisticFamily(name))
}

/**
 * This hook will expose a function for getting and setting whether the
 * statistics are shown. Use this hook if you intend to both read and write the
 * filter value.
 *
 * @param {string} name Name of the filter.
 * @returns Array containing the value and setter function for it.
 */
export function useStatisticState(name) {
  return useRecoilState(statisticFamily(name))
}

/**
 * This hook will expose a function for reading a list of anchored quantities.
 *
 * @returns A list containing the anchored quantity names.
 */
export function useStatisticsValue() {
  return useRecoilValue(statisticsState)
}

export const lockedFamily = atomFamily({
  key: 'lockedFamily',
  default: false
})

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
 * Hook that returns a function for updating the query string.
 *
 * @returns {function} A function that updates the query string to reflect the
 * current search page state.
 */
export function useUpdateQueryString() {
  const history = useHistory()
  const query = useQuery()
  const locked = useRecoilValue(lockedState)
  const statistics = useRecoilValue(statisticsState)

  const updateQueryString = useCallback(() => {
    const queryString = searchToQs(query, locked, statistics)
    history.replace(history.location.pathname + '?' + queryString)
  }, [query, locked, statistics, history])

  return updateQueryString
}

/**
 * Converts a query string into a valid query object.
 *
 * @param {string} queryString URL querystring, encoded or not.
 * @returns Returns an object containing the filters. Values are converted into
 * datatypes that are directly compatible with the filter components.
 */
function qsToSearch(queryString) {
  const queryObj = qs.parse(queryString, {comma: true})

  // Deserialize statistics
  let statistics
  const stats = queryObj.statistics
  if (stats) {
    statistics = {}
    if (isArray(stats)) {
      for (const stat of stats) {
        statistics[stat] = true
      }
    } else {
      statistics[stats] = true
    }
    delete queryObj.statistics
  }

  // Deserialize query
  const query = {}
  for (let [key, value] of Object.entries(queryObj)) {
    const split = key.split(':')
    key = split[0]
    let newKey = filterFullnames[key] || key
    const valueGUI = toGUIFilter(newKey, value)
    if (split.length !== 1) {
      const op = split[1]
      const oldValue = query[newKey]
      if (!oldValue) {
        query[newKey] = {[op]: valueGUI}
      } else {
        query[newKey][op] = valueGUI
      }
    } else {
      query[newKey] = valueGUI
    }
  }
  return [query, statistics]
}

/**
 * Converts a query into an object that uses NOMAD API query string keys and values.
 * @param {object} query A query object representing the currently active
 * filters.
 * @returns {object} An object that can be serialized into a query string
 */
export function searchToQsData(query, locked, statistics) {
  locked = locked || {}
  const queryStringQuery = {}

  // The query is serialized first: locked items will not be displayed in the
  // URL
  if (query) {
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
  }

  // The shown statistics are serialized here: the order is preserved
  if (!isEmpty(statistics)) {
    queryStringQuery.statistics = Object.keys(statistics)
  }

  return queryStringQuery
}

/**
 * Converts a query into a valid query string.
 * @param {object} query A query object representing the currently active
 * filters.
 * @returns URL querystring, not encoded if possible to improve readability.
 */
function searchToQs(query, locked, statistics) {
  const queryData = searchToQsData(query, locked, statistics)
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
export function useAgg(name, update = true, restrict = undefined, delay = 500) {
  const {api} = useApi()
  const { resource } = useSearchContext()
  const [results, setResults] = useState(undefined)
  const initialAggs = useRecoilValue(initialAggsState)
  const finalRestrict = isNil(restrict) ? filterData[name].exclusive : restrict
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
    if (finalRestrict && query && name in query) {
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
  }, [api, name, finalRestrict, resource])

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
  const firstRender = useRef(true)
  const {resource} = useSearchContext()
  const query = useQuery()
  const updateQueryString = useUpdateQueryString()

  const [results, setResults] = useState([])
  const paginationResponse = useRef({})
  const [pagination, setPagination] = useState(initialPagination || {page_size: 10})

  // This callback will call the API with a debounce by the given delay.
  const callAPI = useCallback((query, pagination) => {
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
        updateQueryString()
      })
      .catch(raiseErrors)
  }, [resource, api, raiseErrors, updateQueryString, paginationResponse])

  // Debounced version of callAPI
  const callAPIDebounced = useCallback(debounce(callAPI, delay), [])

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

    if (firstRender.current) {
      callAPI(query, pagination)
      firstRender.current = false
    } else {
      callAPIDebounced(query, pagination)
    }
  }, [callAPI, callAPIDebounced, query, pagination])

  return {
    data: results,
    pagination: combinePagination(pagination, paginationResponse.current),
    setPagination: setPagination
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
function toAPIAgg(aggs, filter, resource, size) {
  const aggSet = filterData[filter].aggSet
  if (aggSet) {
    for (const [key, type] of Object.entries(aggSet)) {
      const name = resource === 'materials' ? materialNames[key.split(':')[0]] : key
      const agg = aggs[name] || {}
      agg[type] = {
        quantity: name,
        size: 200 // Fixed, large value to get a relevant set of results.
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
