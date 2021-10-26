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
import { setToArray, getDatatype, getSerializer, getDeserializer } from '../../utils'
import searchQuantities from '../../searchQuantities'
import { Quantity, getDimension } from '../../units'
import { useErrors } from '../errors'
import { combinePagination } from '../datatable/Datatable'
import InputList from './input/InputList'
import { inputSectionContext } from './input/InputSection'
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
export const labelExperiment = 'Experiment'
export const labelEELS = 'EELS'
export const labelProperties = 'Properties'
export const labelElectronic = 'Electronic'
export const labelVibrational = 'Vibrational'
export const labelMechanical = 'Mechanical'
export const labelSpectroscopy = 'Spectroscopy'
export const labelAuthor = 'Author / Origin'
export const labelAccess = 'Access'
export const labelDataset = 'Dataset'
export const labelIDs = 'IDs'

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
function registerFilterData(name, group, quantity) {
  const agg = quantity.agg
  const value = quantity.value
  const stats = quantity.stats
  const multiple = quantity.multiple === undefined ? true : quantity.multiple
  const exclusive = quantity.exclusive === undefined ? true : quantity.exclusive
  const options = quantity.options
  const unit = quantity.unit
  const dtype = quantity.dtype
  const label = quantity.label
  const description = quantity.description
  const queryMode = quantity.queryMode

  filters.add(name)
  if (group) {
    filterGroups[group]
      ? filterGroups[group].add(name)
      : filterGroups[group] = new Set([name])
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
  data.stats = stats
  data.options = options
  data.unit = unit || searchQuantities[name]?.unit
  data.dtype = dtype || getDatatype(name)
  data.serializerExact = getSerializer(data.dtype, false)
  data.serializerPretty = getSerializer(data.dtype, true)
  data.dimension = getDimension(data.unit)
  data.deserializer = getDeserializer(data.dtype, data.dimension)
  data.label = label
  data.description = description
  data.queryMode = queryMode || 'any'
  filterData[name] = data
}

function registerFilter(name, group, quantity) {
  registerFilterData(name, group, quantity)
}

function registerFilterNested(name, group, quantity, subQuantities) {
  // Register section
  registerFilterData(name, group, quantity)
  filterData[name].nested = true

  // Register section subquantities
  for (let quantity of subQuantities) {
    let subname = `${name}.${quantity.name}`
    registerFilterData(subname, group, quantity)
  }
}

// Presets for different kind of quantities
const termQuantity = {agg: 'terms', stats: listStatConfig}
const noAggQuantity = {stats: listStatConfig}
const nestedQuantity = {}
const noQueryQuantity = {value: {set: () => {}}, multiple: false}
const rangeQuantity = {agg: 'min_max', multiple: false}

// Filters that directly correspond to a metainfo value
registerFilter('results.material.structural_type', labelMaterial, termQuantity)
registerFilter('results.material.functional_type', labelMaterial, termQuantity)
registerFilter('results.material.compound_type', labelMaterial, termQuantity)
registerFilter('results.material.material_name', labelMaterial, termQuantity)
registerFilter('results.material.chemical_formula_hill', labelElements, termQuantity)
registerFilter('results.material.chemical_formula_anonymous', labelElements, termQuantity)
registerFilter('results.material.n_elements', labelElements, {...rangeQuantity, label: 'Number of Elements', queryMode: 'all'})
registerFilter('results.material.symmetry.bravais_lattice', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.crystal_system', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.structure_name', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.strukturbericht_designation', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.space_group_symbol', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.point_group', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.hall_symbol', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.prototype_aflow_id', labelSymmetry, termQuantity)
registerFilter('results.method.method_name', labelMethod, termQuantity)
registerFilter('results.method.simulation.program_name', labelSimulation, termQuantity)
registerFilter('results.method.simulation.program_version', labelSimulation, termQuantity)
registerFilter('results.method.simulation.dft.basis_set_type', labelDFT, termQuantity)
registerFilter('results.method.simulation.dft.core_electron_treatment', labelDFT, termQuantity)
registerFilter('results.method.simulation.dft.xc_functional_type', labelDFT, {...termQuantity, label: 'XC Functional Type'})
registerFilter('results.method.simulation.dft.relativity_method', labelDFT, termQuantity)
registerFilter('results.method.simulation.gw.type', labelGW, {...termQuantity, label: 'GW Type'})
registerFilter('results.method.experiment.eels.detector_type', labelEELS, termQuantity)
registerFilter('results.method.experiment.eels.resolution', labelEELS, rangeQuantity)
registerFilterNested(
  'results.properties.electronic.band_structure_electronic.band_gap',
  labelElectronic,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...rangeQuantity}
  ]
)
registerFilterNested(
  'results.properties.mechanical.bulk_modulus',
  labelMechanical,
  {...nestedQuantity, label: 'Bulk modulus'},
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...rangeQuantity}
  ]
)
registerFilterNested(
  'results.properties.mechanical.shear_modulus',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...rangeQuantity}
  ]
)
registerFilter('external_db', labelAuthor, {...termQuantity, label: 'External Database'})
registerFilter('authors.name', labelAuthor, {...termQuantity, label: 'Author Name'})
registerFilter('upload_create_time', labelAuthor, rangeQuantity)
registerFilter('datasets.dataset_name', labelDataset, {...noAggQuantity, label: 'Dataset Name'})
registerFilter('datasets.doi', labelDataset, {...noAggQuantity, label: 'Dataset DOI'})
registerFilter('entry_id', labelIDs, noAggQuantity)
registerFilter('upload_id', labelIDs, noAggQuantity)
registerFilter('results.material.material_id', labelIDs, noAggQuantity)
registerFilter('datasets.dataset_id', labelIDs, noAggQuantity)
// Visibility: controls the 'owner'-parameter in the API query, not part of the
// query itself.
registerFilter('visibility', labelAccess, noQueryQuantity)
// Restricted: controls whether materials search is done in a restricted mode.
registerFilter('restricted', undefined, noQueryQuantity)
// Exclusive: controls the way elements search is done.
registerFilter('exclusive', undefined, noQueryQuantity)

// In exclusive element query the elements names are sorted and concatenated
// into a single string.
registerFilter(
  'results.material.elements',
  labelElements,
  {
    stats: ptStatConfig,
    agg: 'terms',
    value: {
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
    multiple: true,
    exclusive: false
  }
)
// Electronic properties: subset of results.properties.available_properties
const electronicOptions = {
  band_structure_electronic: {label: 'Band structure'},
  dos_electronic: {label: 'Density of states'}
}
const electronicProps = new Set(Object.keys(electronicOptions))
registerFilter(
  'electronic_properties',
  labelElectronic,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': 'terms'},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => electronicProps.has(value.value)))
    },
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: electronicOptions,
    label: 'Electronic properties',
    description: 'The electronic properties that are present in an entry.'
  }
)
// Vibrational properties: subset of results.properties.available_properties
export const vibrationalOptions = {
  dos_phonon: {label: 'Phonon density of states'},
  band_structure_phonon: {label: 'Phonon band structure'},
  energy_free_helmholtz: {label: 'Helmholtz free energy'},
  heat_capacity_constant_volume: {label: 'Heat capacity constant volume'}
}
const vibrationalProps = new Set(Object.keys(vibrationalOptions))
registerFilter(
  'vibrational_properties',
  labelVibrational,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': 'terms'},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => vibrationalProps.has(value.value)))
    },
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: vibrationalOptions,
    label: 'Vibrational properties',
    description: 'The vibrational properties that are present in an entry.'
  }
)
// Mechanical properties: subset of results.properties.available_properties
export const mechanicalOptions = {
  energy_volume_curve: {label: 'Energy-volume curve'},
  bulk_modulus: {label: 'Bulk modulus'},
  shear_modulus: {label: 'Shear modulus'}
}
const mechanicalProps = new Set(Object.keys(mechanicalOptions))
registerFilter(
  'mechanical_properties',
  labelMechanical,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': 'terms'},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => mechanicalProps.has(value.value)))
    },
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: mechanicalOptions,
    label: 'Mechanical properties',
    description: 'The mechanical properties that are present in an entry.'
  }
)
// Spectroscopic properties: subset of results.properties.available_properties
export const spectroscopicOptions = {
  eels: {label: 'Electron energy loss spectrum'}
}
const spectroscopicProps = new Set(Object.keys(spectroscopicOptions))
registerFilter(
  'spectroscopic_properties',
  labelSpectroscopy,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': 'terms'},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => spectroscopicProps.has(value.value)))
    },
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: spectroscopicOptions,
    label: 'Spectroscopic properties',
    description: 'The spectroscopic properties that are present in an entry.'
  }
)
// EELS energy window: a slider that combines two metainfo values: min_energy
// and max_energy.
registerFilter(
  'results.method.experiment.eels.energy_window',
  labelEELS,
  {
    stats: listStatConfig,
    agg: {
      set: {
        'results.method.experiment.eels.min_energy': 'min_max',
        'results.method.experiment.eels.max_energy': 'min_max'
      },
      get: (aggs) => {
        const min = aggs['results.method.experiment.eels.min_energy']
        const max = aggs['results.method.experiment.eels.max_energy']
        return [min.min_max.data[0], max.min_max.data[1]]
      }
    },
    value: {
      set: (newQuery, oldQuery, value) => {
        newQuery['results.method.experiment.eels.min_energy'] = {gte: value.gte}
        newQuery['results.method.experiment.eels.max_energy'] = {lte: value.lte}
      }
    },
    multiple: false,
    exlusive: false,
    options: vibrationalOptions,
    dtype: 'number',
    unit: 'joule',
    label: 'Energy Window',
    description: 'Defines bounds for the minimum and maximum energies in the spectrum.'
  }
)

// The filter abbreviation mapping has to be done only after all filters have
// been registered.
const abbreviations = {}
const nameAbbreviationPairs = [...filters].map(
  fullname => [fullname, fullname.split('.').pop()])
for (const [fullname, abbreviation] of nameAbbreviationPairs) {
  const old = abbreviations[abbreviation]
  if (old === undefined) {
    abbreviations[abbreviation] = 1
  } else {
    abbreviations[abbreviation] += 1
  }
  filterAbbreviations[fullname] = fullname
  filterFullnames[fullname] = fullname
}
for (const [fullname, abbreviation] of nameAbbreviationPairs) {
  if (abbreviations[abbreviation] === 1) {
    filterAbbreviations[fullname] = abbreviation
    filterFullnames[abbreviation] = fullname
  }
}

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
      setLocked(toGUIFilter(filtersLocked))
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
    useIsCollapsed: () => useRecoilValue(isCollapsedState),
    useSetIsCollapsed: () => useSetRecoilState(isCollapsedState),
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
const isMenuOpenState = atom({
  key: 'isMenuOpen',
  default: false
})
const isCollapsedState = atom({
  key: 'isCollapsed',
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
  // See in which context this filter is being called. If defined within an
  // inputSectionContext, we set the filter inside nested query.
  const sectionContext = useContext(inputSectionContext)
  const section = sectionContext?.section
  const subname = useMemo(() => name.split('.').pop(), [name])

  const value = useRecoilValue(queryFamily(section || name))
  return section
    ? value?.[subname]
    : value
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
  // See in which context this filter is being called. If defined within an
  // inputSectionContext, we set the filter inside nested query.
  const sectionContext = useContext(inputSectionContext)
  const section = sectionContext?.section
  const subname = useMemo(() => name.split('.').pop(), [name])

  const setter = useSetRecoilState(queryFamily(section || name))

  const handleSet = useCallback((value) => {
    section
      ? setter(old => {
        const newValue = isNil(old) ? {} : {...old}
        newValue[subname] = value
        return newValue
      })
      : setter(value)
  }, [section, subname, setter])

  return handleSet
}

/**
 * This hook will expose a function for getting and setting filter values. Use
 * this hook if you intend to both read and write the filter value.
 *
 * @param {string} name Name of the filter.
 * @returns Array containing the filter value and setter function for it.
 */
export function useFilterState(name, section) {
  const value = useFilterValue(name, section)
  const setter = useSetFilter(name, section)
  return useMemo(() => [value, setter], [value, setter])
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
  const query = toGUIFilter(queryObj)
  return [query, statistics]
}

/**
 * Used to create an object that represents the current search context state in
 * a serializable format. Can e.g. be used to build a query string.
 *
 * @param {object} search Object representing the currently active search
 * context.
 *  - query: Object representing the active search filters.
 *  - locked: Object representing the currently locked filters.
 *  - statistics: Object containing the currently shown statistics
 * @returns {object} An object that can e.g. be serialized into a query string.
 */
export function searchToQsData(search) {
  const query = search.query
  const locked = search.locked || {}
  const statistics = search.statistics

  // Used to recursively convert the query into a serializable format.
  function convert(key, value, path) {
    if (locked[key]) {
      return undefined
    }
    // If the key is an operator, the filter name is read from the path.
    const opKeys = new Set(['lte', 'lt', 'gte', 'gt'])
    const fullPath = path ? `${path}.${key}` : key
    const filterPath = opKeys.has(key) ? path : fullPath
    let newValue
    if (isPlainObject(value)) {
      newValue = {}
      for (let [keyInner, valueInner] of Object.entries(value)) {
        const valueConverted = convert(
          keyInner,
          valueInner,
          path ? fullPath : key
        )
        newValue[keyInner] = valueConverted
      }
    } else {
      const serializer = filterData[filterPath].serializerExact
      if (isArray(value)) {
        newValue = value.map(serializer)
      } else if (value instanceof Set) {
        newValue = [...value].map(serializer)
      } else {
        newValue = serializer(value)
      }
    }
    return newValue
  }

  // The query is serialized first: locked items will not be displayed in the
  // URL
  const queryStringQuery = {}
  if (query) {
    for (const [key, value] of Object.entries(query)) {
      const valueConverted = convert(key, value)
      if (!isNil(valueConverted)) {
        // const newKey = filterAbbreviations[key] || key
        const newKey = key
        queryStringQuery[newKey] = valueConverted
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
  const queryData = searchToQsData({query, locked, statistics, abbreviate: true})
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
  const agg = useMemo(() => (
    aggs && {total: aggs.total, data: aggs.data?.[name]}
  ), [aggs, name])
  return agg
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
    // be returned. Notice that the original query is made immutable by
    // Recoil.js: we have to work on a copy.
    let queryCleaned = {...query}
    if (finalRestrict && query) {
      if (name in query) {
        delete queryCleaned[name]
      }
      const splitted = name.split('.')
      const sectionName = splitted.slice(0, -1).join('.')
      const propName = splitted.pop()
      if (sectionName in query && propName in query[sectionName]) {
        const section = {...query[sectionName]}
        delete section[propName]
        queryCleaned[sectionName] = section
      }
    }
    queryCleaned = toAPIFilter(queryCleaned, resource, query.restricted)
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
  const callAPI = useCallback((query, pagination, updateQueryString) => {
    const isExtend = pagination.page_after_value

    const request = {
      owner: query.visibility || 'visible',
      query: toAPIFilter(query, resource, query.restricted),
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
  }, [resource, api, raiseErrors, paginationResponse])

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
      callAPI(query, pagination, updateQueryString)
      firstRender.current = false
    } else {
      callAPIDebounced(query, pagination, updateQueryString)
    }
  }, [callAPI, callAPIDebounced, query, pagination, updateQueryString])

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
export function toAPIFilter(query, resource, restricted) {
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

  // Create the API-compatible keys and values.
  let queryNormalized = {}
  for (const [k, v] of Object.entries(queryCustomized)) {
    const [newKey, newValue] = toAPIFilterSingle(k, v)
    const splitted = newKey.split(':')
    const filterName = splitted[0]
    const queryMode = splitted.length > 1 ? splitted[1] : undefined
    let finalKey = resource === 'materials' ? materialNames[filterName] : filterName
    finalKey = queryMode ? `${finalKey}:${queryMode}` : finalKey
    queryNormalized[finalKey] = newValue
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
 * Cleans a single filter value into a form that is supported by the API. This includes:
 * - Sets are transformed into Arrays
 * - Quantities are converted to SI values.
 *
 * @returns {any} The filter value in a format that is suitable for the API.
 */
function toAPIFilterSingle(key, value, path = undefined) {
  // Determine the API-compatible value.
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
  } else if (isPlainObject(value)) {
    newValue = {}
    for (let [keyInner, valueInner] of Object.entries(value)) {
      const [apiKey, apiValue] = toAPIFilterSingle(keyInner, valueInner, key)
      if (!isNil(apiValue)) {
        newValue[apiKey] = apiValue
      }
    }
    if (isEmpty(newValue)) {
      newValue = undefined
    }
  } else {
    newValue = value
  }

  // Determine the final API key. It depends on the particular queryMode.
  let queryMode
  if (isArray(newValue)) {
    const fullPath = path ? `${path}.${key}` : key
    queryMode = filterData[fullPath]?.queryMode
  }
  const newKey = queryMode ? `${key}:${queryMode}` : key

  return [newKey, newValue]
}

/**
 * Cleans an entire query object into a form that is supported by the GUI. This includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 *
 * @returns {any} The filter object in a format that is suitable for the GUI.
 */
export function toGUIFilter(query, units = undefined) {
  const newQuery = {}
  for (let [key, value] of Object.entries(query)) {
    let newKey = filterFullnames[key] || key
    const valueGUI = toGUIFilterSingle(newKey, value, units)
    newQuery[newKey] = valueGUI
  }
  return newQuery
}

/**
 * Cleans a single filter value into a form that is supported by the GUI. This includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 *
 * @returns {any} The filter value in a format that is suitable for the GUI.
 */
export function toGUIFilterSingle(key, value, units = undefined, path = undefined) {
  let newValue
  const fullPath = path ? `${path}.${key}` : key
  if (isPlainObject(value)) {
    newValue = {}
    for (let [keyInner, valueInner] of Object.entries(value)) {
      const valueConverted = toGUIFilterSingle(keyInner, valueInner, units, fullPath)
      if (!isNil(valueConverted)) {
        newValue[keyInner] = valueConverted
      }
    }
  } else {
    // If the key is an operator, the filter name is read from the path.
    const opKeys = new Set(['lte', 'lt', 'gte', 'gt'])
    const filterPath = opKeys.has(key) ? path : fullPath
    let multiple = filterData[filterPath].multiple
    const deserializer = filterData[filterPath].deserializer
    if (isArray(value)) {
      newValue = new Set(value.map((v) => deserializer(v, units)))
    } else {
      newValue = deserializer(value, units)
    }
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
