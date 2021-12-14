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
  isSet
} from 'lodash'
import qs from 'qs'
import PropTypes from 'prop-types'
import { useHistory } from 'react-router-dom'
import { useApi } from '../api'
import { setToArray } from '../../utils'
import { Quantity } from '../../units'
import { useErrors } from '../errors'
import { combinePagination } from '../datatable/Datatable'
import { inputSectionContext } from './input/InputSection'
import {
  filterDataGlobal,
  filterAbbreviations,
  filterFullnames,
  materialNames,
  entryNames
} from './FilterRegistry'

/**
 * React context that provides access to the search state implemented with
 * Recoil.js. The purpose of this context is to hide the Recoil.js
 * implementation details and provide a clean access to individual states.
 *
 * Each search filter, aggregation query and aggregation result is mapped into a
 * separate Recoil.js Atom. This allows components to hook into individual
 * search parameters (both for setting and reading their value). This performs
 * much better than having larger Atoms, as this would cause several components
 * to render even if they are not affected by some other search filter or
 * result.
 *
 * Access to the search state is typically made through functions that can be
 * accessed from the context. This way components can hook into individual
 * values of information and do not have to re-render. If we were to provide
 * e.g. the filter values directly within the search context, each component
 * depending on the context would re-render for each filter change regardless if
 * it actually changes the state of that component or not.
 */
const orderByMap = {
  'entries': 'upload_create_time',
  'materials': 'chemical_formula_hill'
}
let indexContext = 0
let indexFilters = 0
let indexLocked = 0

export const searchContext = React.createContext()
export const SearchContext = React.memo(({
  resource,
  filtersLocked,
  initialPagination,
  children
}) => {
  const {api} = useApi()
  const {raiseErrors} = useErrors()
  const oldQuery = useRef(undefined)
  const oldPagination = useRef(undefined)
  const paginationResponse = useRef(undefined)
  const updatedAggsMap = useRef({})
  const apiQueue = useRef([])
  const apiMap = useRef({})
  const updatedFilters = useRef(new Set())
  const refreshFilters = useRef(new Set())
  const firstLoad = useRef(true)

  // Initialize the set of available filters. This may depend on the resource.
  const [filtersLocal, filterDataLocal] = useMemo(() => {
    const filtersLocal = new Set()
    const filterDataLocal = []
    for (const [key, value] of Object.entries(filterDataGlobal)) {
      if (value.resources.has(resource)) {
        filtersLocal.add(key)
        filterDataLocal[key] = value
      }
    }
    return [filtersLocal, filterDataLocal]
  }, [resource])
  const filters = useState(filtersLocal)[0]
  const filterData = useState(filterDataLocal)[0]

  // Initialize the search context state: any parameters in the URL are read and
  // default values as specified in filter registry are loaded
  const [initialQuery, initialAggs, initialStatistics, filterDefaults] = useMemo(() => {
    const filterDefaults = {}
    for (const [key, value] of Object.entries(filterData)) {
      if (!isNil(value.default)) {
        filterDefaults[key] = value.default
      }
    }
    const location = window.location.href
    const split = location.split('?')
    let queryURL = {}
    let statisticsURL = {}
    if (split.length !== 1) {
      const qs = split.pop();
      [queryURL, statisticsURL] = qsToSearch(qs)
    }

    const initialAggs = {}
    for (const key of filters) {
      initialAggs[key] = {
        default: {update: false}
      }
    }
    return [
      {...queryURL, ...toGUIFilter(filtersLocked)},
      initialAggs,
      statisticsURL,
      filterDefaults
    ]
  }, [filterData, filters, filtersLocked])

  // Initialize a bunch of Recoil.js states and hooks. Notice how we are not using a set
  // of global states, but instead each SearchContext gets it's own states. This
  // way the contexts stay separated and several of them can be active at once.
  const [
    useFilterLocked,
    useFiltersLocked,
    useFiltersLockedState,
    useFilterValue,
    useSetFilter,
    useFilterState,
    useFiltersState,
    useResetFilters,
    useUpdateQueryString,
    queryState,
    aggsState,
    paginationState,
    resultsState,
    aggsResponseState,
    isMenuOpenState,
    isCollapsedState,
    isStatisticsEnabledState,
    useStatisticValue,
    useSetStatistic,
    useStatisticState,
    useStatisticsValue,
    useResults,
    useAgg,
    useSetFilters
  ] = useMemo(() => {
    const queryFamily = atomFamily({
      key: `queryFamily_${indexContext}`,
      default: (name) => initialQuery[name]
    })
    // Used to get/set the locked state of all filters at once
    const filtersState = selector({
      key: `filtersState_${indexContext}`,
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

    const guiLocked = toGUIFilter(filtersLocked)
    const lockedFamily = atomFamily({
      key: `lockedFamily_${indexContext}`,
      default: (name) => !isNil(guiLocked?.[name])
    })

    // Used to set the locked state of several filters at once
    const lockedState = selector({
      key: `lockedState_${indexContext}`,
      get: ({get}) => {
        const locks = {}
        for (let key of filters) {
          const filter = get(lockedFamily(key))
          locks[key] = filter
        }
        return locks
      }
    })

    /**
     * A Recoil.js selector that aggregates all the currently set filters into a
     * single query object used by the API.
     */
    const queryState = selector({
      key: `query_${indexContext}`,
      get: ({get}) => {
        let query = {}
        for (let key of filters) {
          const filter = get(queryFamily(key))
          if (filter !== undefined) {
            query[key] = filter
          }
        }
        return query
      },
      set: ({ set }, data) => {
        for (let filter of filters) {
          set(queryFamily(filter), undefined)
        }
        if (data) {
          for (const [key, value] of Object.entries(data)) {
            set(queryFamily(key), value)
          }
        }
      }
    })

    const isStatisticsEnabledState = atom({
      key: `statisticsEnabled_${indexContext}`,
      default: true
    })
    const isMenuOpenState = atom({
      key: `isMenuOpen_${indexContext}`,
      default: false
    })
    const isCollapsedState = atom({
      key: `isCollapsed_${indexContext}`,
      default: false
    })

    const paginationState = atom({
      key: `pagination_${indexContext}`,
      default: initialPagination || {
        order_by: orderByMap[resource],
        page_size: 20
      }
    })

    const statisticFamily = atomFamily({
      key: `statisticFamily_${indexContext}`,
      default: (name) => initialStatistics[name]
    })

    // Used to get/set the the statistics configuration of all filters
    const statisticsState = selector({
      key: `statisticsState_${indexContext}`,
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
    function useStatisticValue(name) {
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
    function useSetStatistic(name) {
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
    function useStatisticState(name) {
      return useRecoilState(statisticFamily(name))
    }

    /**
     * This hook will expose a function for reading a list of anchored quantities.
     *
     * @returns A list containing the anchored quantity names.
     */
    function useStatisticsValue() {
      return useRecoilValue(statisticsState)
    }

    const resultsState = atom({
      key: `results_${indexContext}`,
      default: {
        pagination: {}
      }
    })

    const aggsFamily = atomFamily({
      key: `aggsFamily_${indexContext}`,
      default: (name) => initialAggs[name]
    })

    /**
     * A Recoil.js selector that aggregates all the currently set aggregation
     * requests into a single query object used by the API.
     */
    const aggsState = selector({
      key: `aggs_${indexContext}`,
      get: ({get}) => {
        let query = {}
        for (let key of filters) {
          const filter = get(aggsFamily(key))
          if (filter !== undefined) {
            query[key] = filter
          }
        }
        return query
      }
    })

    // Atom for each aggregation response.
    const aggsResponseFamily = atomFamily({
      key: `aggsResponseFamily_${indexContext}`,
      default: undefined
    })

    // Recoil.js selector for updating the aggs response in one go.
    const aggsResponseState = selector({
      key: `aggsResponse_${indexContext}`,
      set: ({ set }, data) => {
        if (data) {
          for (const [key, value] of Object.entries(data)) {
            set(aggsResponseFamily(key), value)
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
    const useFilterValue = (name) => {
      // See in which context this filter is being called. If defined within an
      // inputSectionContext, we set the filter inside nested query.
      const sectionContext = useContext(inputSectionContext)
      const section = sectionContext?.section
      const nested = sectionContext?.nested
      const subname = useMemo(() => name.split('.').pop(), [name])

      const value = useRecoilValue(queryFamily(section || name))
      return (section && nested)
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
    const useSetFilter = (name) => {
      // See in which context this filter is being called. If defined within an
      // inputSectionContext, we set the filter inside nested query.
      const sectionContext = useContext(inputSectionContext)
      const section = sectionContext?.section
      const nested = sectionContext?.nested
      const subname = useMemo(() => name.split('.').pop(), [name])

      const setter = useSetRecoilState(queryFamily(section || name))

      const handleSet = useCallback((value) => {
        updatedFilters.current.add(name)
        section && nested
          ? setter(old => {
            const newValue = isNil(old) ? {} : {...old}
            newValue[subname] = value
            return newValue
          })
          : setter(value)
      }, [section, subname, setter, name, nested])

      return handleSet
    }

    /**
     * This hook will expose a function for getting and setting filter values. Use
     * this hook if you intend to both read and write the filter value.
     *
     * @param {string} name Name of the filter.
     * @returns Array containing the filter value and setter function for it.
     */
    const useFilterState = (name, section) => {
      const value = useFilterValue(name, section)
      const setter = useSetFilter(name, section)
      return useMemo(() => [value, setter], [value, setter])
    }

    /**
     * Hook that returns a function for updating the query string.
     *
     * @returns {function} A function that updates the query string to reflect the
     * current search page state.
     */
    const useUpdateQueryString = () => {
      const history = useHistory()
      const query = useRecoilValue(queryState)
      const locked = useRecoilValue(lockedState)
      const statistics = useRecoilValue(statisticsState)

      return useCallback(() => {
        const queryString = searchToQs(query, locked, statistics)
        history.replace(history.location.pathname + '?' + queryString)
      }, [history, locked, query, statistics])
    }

    /**
     * Returns a function that can be called to reset all current filters.
     *
     * @returns Function for resetting all filters.
     */
    const useResetFilters = () => {
      const locked = useRecoilValue(lockedState)
      const reset = useRecoilCallback(({set}) => () => {
        for (let filter of filters) {
          if (!locked[filter]) {
            set(queryFamily(filter), undefined)
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
    const useFilterLocked = (name) => {
      return useRecoilValue(lockedFamily(name))
    }

    /**
     * This hook will expose a function for reading the locked status of all
     * filters.
     *
     * @returns An object containing a mapping from filter name to a boolean
     * indicating whether it is locked or not.
     */
    const useFiltersLocked = () => {
      return useRecoilValue(lockedState)
    }

    /**
     * This hook will expose a function for reading if the given set of filters are
     * locked.
     *
     * @param {string} names Names of the filters.
     * @returns Array containing the filter values in a map and a setter function.
     */
    const useFiltersLockedState = (names) => {
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

    /**
     * This hook will expose a function for getting and setting filter values for
     * the specified list of filters. Use this hook if you intend to both read and
     * write the filter values.
     *
     * @param {string} names Names of the filters.
     * @returns Array containing the filter values in a map and a setter function.
     */
    const useFiltersState = (names) => {
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
     * Hook for returning an object containing the currently fetched results,
     * pagination and a function for changing the pagination.
     *
     * @returns {object} {data, pagination, setPagination}
     */
    const useResults = () => useRecoilValue(resultsState)
    /**
     * Hook for modifying an aggregation request and fetching the latest values for
     * this aggregation.
     *
     * @param {string} name The filter name
     * @param {bool} update Whether the hook needs to react to changes in the
     * current query context. E.g. if the component showing the data is not visible,
     * this can be set to false.
     *
     * @returns {object} An object containing the aggregation results: the layout is
     * specific for each aggregation type.
     */

    const useAgg = (name, update = true, size = undefined, id = 'default') => {
      const setAgg = useSetRecoilState(aggsFamily(name))
      const aggResponse = useRecoilValue(aggsResponseFamily(name))
      const aggSize = size || filterData[name].aggSize

      useEffect(() => {
        setAgg(old => {
          const config = {update, size: aggSize}
          const newAgg = old ? {...old, [id]: config} : {[id]: config}
          return newAgg
        })
      }, [name, update, aggSize, id, setAgg])

      return aggResponse
    }

    /**
     * This hook will expose a function for setting the values of all filters.
     *
     * @returns An object containing a mapping from filter name to a boolean
     * indicating whether it is locked or not.
     */
    const useSetFilters = () => {
      return useSetRecoilState(filtersState)
    }

    ++indexContext
    return [
      useFilterLocked,
      useFiltersLocked,
      useFiltersLockedState,
      useFilterValue,
      useSetFilter,
      useFilterState,
      useFiltersState,
      useResetFilters,
      useUpdateQueryString,
      queryState,
      aggsState,
      paginationState,
      resultsState,
      aggsResponseState,
      isMenuOpenState,
      isCollapsedState,
      isStatisticsEnabledState,
      useStatisticValue,
      useSetStatistic,
      useStatisticState,
      useStatisticsValue,
      useResults,
      useAgg,
      useSetFilters
    ]
  }, [initialQuery, filters, filtersLocked, initialStatistics, initialAggs, initialPagination, resource, filterData])

  const setResults = useSetRecoilState(resultsState)
  const updateAggsResponse = useSetRecoilState(aggsResponseState)
  const aggs = useRecoilValue(aggsState)
  const query = useRecoilValue(queryState)
  const [pagination, setPagination] = useRecoilState(paginationState)
  const updateQueryString = useUpdateQueryString()
  const isStatisticsEnabled = useRecoilValue(isStatisticsEnabledState)

  // All of the heavier pre-processing, checking, etc. should be done in this
  // function, as it is the final one that gets called after the debounce
  // interval.
  const apiCall = useCallback((query, aggs, pagination, queryChanged, paginationChanged, updateAggs, aggNames, refresh = false) => {
    // Create the final search object.
    let apiQuery = {...query}
    if (filterDefaults) {
      for (const [key, value] of Object.entries(filterDefaults)) {
        if (isNil(query[key])) {
          apiQuery[key] = value
        }
      }
    }
    const search = {
      owner: apiQuery.visibility,
      query: toAPIFilter(apiQuery, resource, filterDefaults),
      aggregations: toAPIAgg(
        aggs,
        aggNames,
        refresh ? refreshFilters.current : updatedFilters.current,
        resource
      ),
      pagination: {...pagination}
    }
    // When aggregations have changed but the query has not, we request only the
    // aggregation data without any hits.
    if (updateAggs && !queryChanged && !paginationChanged) {
      search.pagination = {page_size: 0}
      search.required = { include: [] }
    }

    // If query changes or the pagination changes without page_after_value
    // changing, we request the first page. This is a simple way ensure that
    // results are consistent.
    if (queryChanged || (paginationChanged && (pagination.page_after_value && pagination.page_after_value === paginationResponse.current.page_after_value))) {
      search.pagination.page_after_value = undefined
      search.pagination.next_page_after_value = undefined
    }

    // Due to the queueing mechanism we can now already update the reference to
    // contain the latest information about what was updated by this query. This
    // ensures that the next query immediately knows the current state even
    // before the API call is finished. When query changes, the aggregations are
    // reset to only contain the ones fetched simultaneously with the query.
    // Otherwise any new aggregations are just added to the set of already
    // updated ones. The list of updated filters are always reset, and the list
    // of updated filters is stored for refresh purposes.
    if (queryChanged) {
      const mapping = {}
      aggNames.forEach((agg) => { mapping[agg] = aggs[agg] })
      updatedAggsMap.current = mapping
    } else {
      aggNames.forEach((agg) => { updatedAggsMap.current[agg] = aggs[agg] })
    }
    refreshFilters.current = updatedFilters.current
    updatedFilters.current = new Set()
    firstLoad.current = false
    oldQuery.current = query
    oldPagination.current = pagination

    // As we cannot guarantee the order in which the API calls finish, we push
    // all calls into a queue. The API calls are always made instantly, but the
    // queue makes sure that API calls get resolved in the original order not
    // matter how long the actual call takes.
    function resolve(prop) {
      const {response, timestamp, queryChanged, paginationChanged, aggNames, search, resource} = prop
      let next = apiQueue.current[0]
      if (next !== timestamp) {
        apiMap.current[timestamp] = prop
        return
      }
      // Update the aggregations if new aggregation data is received. The old
      // aggregation data is preserved and new information is updated.
      if (!isEmpty(response.aggregations)) {
        const newAggs = toGUIAgg(response.aggregations, aggNames, resource)
        updateAggsResponse(newAggs)
      }
      // Update the query results if new data is received.
      if (queryChanged || paginationChanged) {
        const isExtend = search.pagination.page_after_value
        paginationResponse.current = response.pagination
        setResults(old => {
          const newResults = old ? {...old} : {}
          isExtend ? newResults.data = [...newResults.data, ...response.data] : newResults.data = response.data
          newResults.pagination = combinePagination(search.pagination, response.pagination)
          newResults.setPagination = setPagination
          return newResults
        })
      }
      // Remove this query from queue and see if next can be resolved.
      apiQueue.current.shift()
      const nextTimestamp = apiQueue.current[0]
      const nextResolve = apiMap.current[nextTimestamp]
      if (nextResolve) {
        resolve(nextResolve)
      }
    }
    const timestamp = Date.now()
    apiQueue.current.push(timestamp)
    api.query(resource, search, true).then((response) => resolve({
      response,
      timestamp,
      queryChanged,
      paginationChanged,
      aggNames,
      search,
      resource
    })).catch(raiseErrors)
  }, [filterDefaults, resource, api, raiseErrors, updateAggsResponse, setResults, setPagination])

  // This is a debounced version of apiCall.
  const apiCallDebounced = useCallback(debounce(apiCall, 400), [])

  // When query, aggregation or pagination changes, make an API call. The API
  // call is made immediately on first render. On subsequent renders it will be
  // debounced.
  useEffect(() => {
    // If the query and pagination has not changed AND aggregations do not need
    // to be updated, no update is necessary.
    const queryChanged = query !== oldQuery.current
    const paginationChanged = pagination !== oldPagination.current
    let aggNames
    const reducedAggs = reduceAggs(aggs, updatedAggsMap.current, queryChanged)
    if (!isStatisticsEnabled) {
      aggNames = []
    } else {
      aggNames = Object.keys(reducedAggs).filter((key) => reducedAggs[key].update)
    }
    const updateAggs = aggNames.length > 0
    if (!paginationChanged && !queryChanged && !updateAggs) {
      return
    }

    // The API calls is made immediately when requesting the first set of
    // results, when the pagination changes or when only aggregations need to be
    // updated. Otherwise it is debounced.
    if (firstLoad.current || paginationChanged || !queryChanged) {
      apiCall(query, reducedAggs, pagination, queryChanged, paginationChanged, updateAggs, aggNames)
    } else {
      apiCallDebounced(query, reducedAggs, pagination, queryChanged, paginationChanged, updateAggs, aggNames)
    }
  }, [query, aggs, pagination, apiCall, apiCallDebounced, isStatisticsEnabled])

  // Hook for refreshing the results
  const useRefresh = useCallback(() => {
    const query = useRecoilValue(queryState)
    const aggs = useRecoilValue(aggsState)
    const pagination = useRecoilValue(paginationState)
    const queryChanged = true
    const paginationChanged = false
    const updateAggs = true
    const aggNames = Object.keys(aggs)

    const refresh = useCallback(() => {
      apiCallDebounced(query, aggs, pagination, queryChanged, paginationChanged, updateAggs, aggNames, true)
    }, [aggNames, aggs, pagination, paginationChanged, query, queryChanged, updateAggs])
    return refresh
  }, [aggsState, apiCallDebounced, paginationState, queryState])

  // This updated the query string to represent the latest value within the
  // search context.
  useEffect(() => {
    updateQueryString()
  }, [updateQueryString])

  // The context contains a set of functions that can be used to hook into
  // different data.
  const values = useMemo(() => ({
    resource,
    useIsMenuOpen: () => useRecoilValue(isMenuOpenState),
    useSetIsMenuOpen: () => useSetRecoilState(isMenuOpenState),
    useIsCollapsed: () => useRecoilValue(isCollapsedState),
    useSetIsCollapsed: () => useSetRecoilState(isCollapsedState),
    useIsStatisticsEnabled: () => useRecoilValue(isStatisticsEnabledState),
    useSetIsStatisticsEnabled: () => useSetRecoilState(isStatisticsEnabledState),
    useQuery: () => useRecoilValue(queryState),
    useFilterValue: useFilterValue,
    useSetFilter: useSetFilter,
    useFilterState: useFilterState,
    useFiltersState: useFiltersState,
    useResetFilters: useResetFilters,
    useRefresh: useRefresh,
    useFilterLocked: useFilterLocked,
    useFiltersLocked: useFiltersLocked,
    useFiltersLockedState: useFiltersLockedState,
    useStatisticValue: useStatisticValue,
    useSetStatistic: useSetStatistic,
    useStatisticState: useStatisticState,
    useStatisticsValue: useStatisticsValue,
    useUpdateQueryString: useUpdateQueryString,
    useResults: useResults,
    useAgg: useAgg,
    useSetFilters: useSetFilters,
    filters: filters,
    filterData: filterData
  }), [
    resource,
    useFilterValue,
    useSetFilter,
    useFilterState,
    useFiltersState,
    useResetFilters,
    useFilterLocked,
    useFiltersLocked,
    useFiltersLockedState,
    useStatisticValue,
    useSetStatistic,
    useStatisticState,
    useStatisticsValue,
    useUpdateQueryString,
    useRefresh,
    useResults,
    useAgg,
    useSetFilters,
    isMenuOpenState,
    isCollapsedState,
    isStatisticsEnabledState,
    queryState,
    filters,
    filterData
  ])

  return <searchContext.Provider value={values}>
    {children}
  </searchContext.Provider>
})
SearchContext.propTypes = {
  resource: PropTypes.string,
  filtersLocked: PropTypes.object,
  initialPagination: PropTypes.object,
  children: PropTypes.node
}

/**
 * Hook for accessing the current SearchContext.
 */
export function useSearchContext() {
  return useContext(searchContext)
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
  let statistics = {}
  const stats = queryObj.statistics
  if (stats) {
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
    const filterData = filterDataGlobal[filterPath]
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
      const serializer = filterData.serializerExact
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
        const newKey = filterAbbreviations[key] || key
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
 * @param {object} query Query object representing the currently active
 * @param {object} locked Object containing the locked status of quantities.
 * @param {object} statistics Object containing which filter statistics are docked.
 * filters.
 * @returns URL querystring, not encoded if possible to improve readability.
 */
function searchToQs(query, locked, statistics) {
  const queryData = searchToQsData({query, locked, statistics, abbreviate: true})
  return qs.stringify(queryData, {indices: false, encode: false})
}

/**
 * Converts the contents of a query into a format that is suitable for the API.
 *
 * Should only be called when making the final API call, as during the
 * construction of the query it is much more convenient to store filters within
 * e.g. Sets.
 *
 * @param {number} query The query object.
 * @param {string} resource The resource we are looking at: entries or materials.
 *
 * @returns {object} A copy of the object with certain items cleaned into a
 * format that is supported by the API.
 */
export function toAPIFilter(query, resource, filterDefaults) {
  let queryCustomized = {}
  if (!query) {
    return undefined
  }

  // Perform custom transformations
  const combine = query.combine
  for (let [k, v] of Object.entries(query)) {
    const data = filterDataGlobal[k]
    const guiOnly = data?.guiOnly
    const setter = data?.valueSet
    if (guiOnly) {
      continue
    }
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
    // When combining results, we split each filter and each filter value into
    // it's own separate entries query. These queries are then joined with
    // 'and'.
    if (combine) {
      const entrySearch = []
      for (const [k, v] of Object.entries(queryNormalized)) {
        if (k.startsWith('entries.')) {
          const splitted = k.split(':')
          const [newKey, queryMode] = splitted.length === 2 ? splitted : [splitted[0], undefined]
          // When the queryMode is 'all', each value can come from a separate
          // entry.
          if (isArray(v) && queryMode === 'all') {
            for (const item of v) {
              entrySearch.push({[newKey]: item})
            }
          } else {
            entrySearch.push({[k]: v})
          }
          delete queryNormalized[k]
        }
      }
      if (entrySearch.length > 0) {
        queryNormalized.and = entrySearch
      }
    // When combining results is not allowed, we simply make a nested query by
    // moving all method/properties filters inside a single entries-subsection.
    } else {
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
    }
  }

  return queryNormalized
}

/**
 * Cleans a single filter value into a form that is supported by the API. This includes:
 * - Sets are transformed into Arrays
 * - Quantities are converted to SI values.
 *
 * @param {string} key Filter name
 * @param {any} value Filter value
 * @param {string} path The full path of the filter.
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
    queryMode = filterDataGlobal[fullPath]?.queryMode
  }
  const newKey = queryMode ? `${key}:${queryMode}` : key

  return [newKey, newValue]
}

/**
 * Cleans an entire query object into a form that is supported by the GUI. This
 * includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 *
 * @param {object} query Query object to transform.
 * @param {object} units The desired unit system used when reading quantities.
 *
 * @returns {any} The filter object in a format that is suitable for the GUI.
 */
export function toGUIFilter(query, units = undefined) {
  const newQuery = {}
  if (query) {
    for (let [key, value] of Object.entries(query)) {
      let newKey = filterFullnames[key] || key
      const valueGUI = toGUIFilterSingle(newKey, value, units)
      newQuery[newKey] = valueGUI
    }
  }
  return newQuery
}

/**
 * Cleans a single filter value into a form that is supported by the GUI. This includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 *
 * @param {string} key Name of the filter or an operator name.
 * @param {any} value Value of the filter or operator
 * @param {object} units Unit system for unit conversion
 * @param {string} path The full path for the filter without an operator name.
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
    let multiple = filterDataGlobal[filterPath].multiple
    const deserializer = filterDataGlobal[filterPath].deserializer
    if (isArray(value) || isSet(value)) {
      newValue = new Set(value.map((v) => deserializer(v, units)))
    } else {
      newValue = deserializer(value, units)
      if (multiple) {
        newValue = new Set([newValue])
      }
    }
  }
  return newValue
}

/**
 * Used to transform a GUI aggregation query into a form that is usable by the
 * API.
 *
 * @param {object} aggs The aggregation data as constructed by the GUI.
 * @param {array} aggNames The aggregation names to update
 * @param {string} resource The resource we are looking at: entries or materials.
 * @param {bool} update Whether to force the update of aggregations, overriding
 * the update-attribute of each aggregation.
 *
 * @returns {object} Aggregation query that is usable by the API.
 */
function toAPIAgg(aggs, aggNames, updatedFilters, resource) {
  const apiAggs = {}
  for (const aggName of aggNames) {
    const agg = aggs[aggName]
    const aggSet = filterDataGlobal[aggName].aggSet
    if (aggSet) {
      for (const [key, type] of Object.entries(aggSet)) {
        // If filter has been updated and the filter values are exclusive, the
        // filter is excluded from the aggregation.
        const exclude = updatedFilters.has(key) && filterDataGlobal[key].exclusive
        const name = resource === 'materials' ? materialNames[key.split(':')[0]] : key
        const apiAgg = apiAggs[name] || {}
        apiAgg[type] = {
          quantity: name,
          exclude_from_search: exclude,
          size: agg.size
        }
        apiAggs[name] = apiAgg
      }
    }
  }
  return apiAggs
}

/**
 * Used to transform an API aggregation result into a form that is usable by the
 * GUI.
 *
 * @param {object} aggs The aggregation data as returned by the API.
 * @param {array} filters The set of targeted filters. Needed because the keys
 * in the aggs dictionary may be different due to custom aggregation set/get.
 * @param {string} resource The resource we are looking at: entries or materials.
 *
 * @returns {object} Aggregation result that is usable by the GUI.
 */
function toGUIAgg(aggs, filters, resource) {
  if (isEmpty(aggs)) {
    return {}
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
    const aggGet = filterDataGlobal[name]?.aggGet
    if (aggGet) {
      const agg = aggGet(aggsNormalized)
      aggsCustomized[name] = {
        data: agg,
        // TODO: Could this total be given by the API directly?
        total: agg[0]?.count && agg.reduce((a, b) => a + b.count, 0)
      }
    }
  }
  return aggsCustomized
}

/**
 * Reduces the current agggregation setup into simpler form that contains only
 * two variables: whether to update the aggregation and with what size if one
 * has been specified.
 *
 * @param {object} aggs The current aggregation configuration.
 * @param {object} oldAggs Reduced aggregation config from latest finished query.
 * @param {bool} queryChanged Whether the query has changed.
 *
 * @returns {object} Reduced aggregation config.
 */
function reduceAggs(aggs, oldAggs, queryChanged) {
  const reducedAggs = {}
  for (let [key, agg] of Object.entries(aggs)) {
    let update = false
    let size = 0

    // Loop through the different configs and see if any of them need to be
    // updated and what is the maximum requested agg size.
    for (let config of Object.values(agg)) {
      if (config.update) {
        update = true
      }
      if (!isNil(config.size) && config.size > size) {
        size = config.size
      }
    }

    // Some aggregations require us to load more data than what we are currently
    // showing (e.g. properties lists).
    const sizeOverride = filterDataGlobal[key].aggSizeOverride
    if (sizeOverride) {
      size = sizeOverride
    }

    // If the query has not changed, see if there is an old aggregation which
    // has a size that is at least as big as the currently requested size.
    if (!queryChanged) {
      const oldAgg = oldAggs[key]
      if (oldAgg) {
        if (!isNil(oldAgg.size) && size > oldAgg.size) {
          update = true
        } else {
          update = false
        }
      }
    }
    const newAgg = {update: update}
    if (size) {
      newAgg.size = size
    }
    reducedAggs[key] = newAgg
  }
  return reducedAggs
}
