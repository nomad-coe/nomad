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
  selectorFamily,
  useSetRecoilState,
  useRecoilValue,
  useRecoilState,
  useRecoilCallback
} from 'recoil'
import {
  debounce,
  isEmpty,
  isArray,
  isBoolean,
  isPlainObject,
  isNil,
  isSet,
  isFunction,
  size,
  isEqual
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
  filterData as filterDataGlobal,
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
  'entries': {order_by: 'upload_create_time', order: 'desc'},
  'materials': {order_by: 'chemical_formula_hill', order: 'asc'}
}
let indexContext = 0
let indexFilters = 0
let indexLocked = 0

/**
 * Used to turn empty containers into undefined which is used to indicate that a
 * search filter has not been specified.
 */
function clearEmpty(value) {
  function isEmpty(value) {
    if (isPlainObject(value)) {
      return Object.values(value).every(isEmpty)
    } else if ((isSet(value) || isArray(value)) && size(value) === 0) {
      return true
    }
    return false
  }
  return isEmpty(value) ? undefined : value
}

export const searchContext = React.createContext()
export const SearchContext = React.memo(({
  resource,
  filtersLocked,
  initialPagination,
  children
}) => {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const oldQuery = useRef(undefined)
  const oldPagination = useRef(undefined)
  const paginationResponse = useRef(undefined)
  const updatedAggsMap = useRef({})
  const apiQueue = useRef([])
  const apiMap = useRef({})
  const updatedFilters = useRef(new Set())
  const firstLoad = useRef(true)
  const disableUpdate = useRef(false)

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
    aggsFamily,
    aggsState,
    paginationState,
    resultsUsedState,
    resultsState,
    apiDataState,
    aggsResponseState,
    isMenuOpenState,
    isCollapsedState,
    isStatisticsEnabledState,
    useStatisticValue,
    useSetStatistic,
    useStatisticState,
    useStatisticsValue,
    useStatisticsState,
    useResults,
    useApiData,
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
        for (const key of filters) {
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
        for (const key of filters) {
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
        const query = {}
        for (const key of filters) {
          const filter = get(queryFamily(key))
          if (filter !== undefined) {
            query[key] = filter
          }
        }
        return query
      },
      set: ({ set }, data) => {
        for (const filter of filters) {
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
        page_size: 20,
        ...orderByMap[resource]
      }
    })

    const resultsUsedState = atom({
      key: `resultsUsed_${indexContext}`,
      default: false
    })

    const statisticFamily = atomFamily({
      key: `statisticFamily_${indexContext}`,
      default: (name) => initialStatistics[name]
    })

    // Used to get/set the the statistics configuration of all filters
    const statisticsState = selector({
      key: `statisticsState_${indexContext}`,
      set: ({set}, stats) => {
        stats && Object.entries(stats).forEach(([key, value]) => set(statisticFamily(key), value))
      },
      get: ({get}) => {
        const stats = {}
        for (const filter of filters) {
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
     * @returns Function for setting the value
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
     * @returns Array containing the value and a function for setting it.
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

    /**
     * This hook will expose a function for reading and writing the object
     * containing all the current statistics.
     *
     * @returns Array containing the value and a function for setting it.
     */
    function useStatisticsState() {
      return useRecoilState(statisticsState)
    }

    const resultsState = atom({
      key: `results_${indexContext}`,
      default: {
        pagination: {}
      }
    })

    const apiDataState = atom({
      key: `apiData_${indexContext}`,
      default: null
    })

    const aggsFamilyRaw = atomFamily({
      key: `aggsFamilyRaw_${indexContext}`,
      default: (name) => initialAggs[name]
    })

    const aggKeys = atom({
      key: `aggKeys_${indexContext}`,
      default: []
    })

    const aggsFamily = selectorFamily({
      key: `aggsFamily_${indexContext}`,
      get: (id) => ({ get }) => {
        return get(aggsFamilyRaw(id))
      },
      set: (id) => ({set}, meal) => {
        set(aggsFamilyRaw(id), meal)
        set(aggKeys, prev => [...prev, id])
      }
    })

    /**
     * A Recoil.js selector that aggregates all the currently set aggregation
     * requests into a single query object used by the API.
     */
    const aggsState = selector({
      key: `aggs_${indexContext}`,
      get: ({get}) => {
        const aggs = {}
        for (const key of get(aggKeys)) {
          const agg = get(aggsFamily(key))
          if (agg !== undefined) {
            aggs[key] = agg
          }
        }
        return aggs
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
      const subname = useMemo(() => section ? name.slice(section.length + 1) : undefined, [name, section])

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
    const useSetFilter = (name) => {
      // See in which context this filter is being called. If defined within an
      // inputSectionContext, we set the filter inside nested query.
      const sectionContext = useContext(inputSectionContext)
      const section = sectionContext?.section
      const subname = useMemo(() => section ? name.slice(section.length + 1) : undefined, [name, section])

      const setter = useSetRecoilState(queryFamily(section || name))

      const handleSet = useCallback((value) => {
        updatedFilters.current.add(name)
        section
          ? setter(old => {
            const newValue = isNil(old) ? {} : {...old}
            if (isFunction(value)) {
              value = value(newValue[subname])
            }
            newValue[subname] = value
            return clearEmpty(newValue)
          })
          : setter(isFunction(value)
            ? (old) => clearEmpty(value(old))
            : clearEmpty(value))
      }, [section, subname, setter, name])
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
        for (const filter of filters) {
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
            for (const key of names) {
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
            for (const key of names) {
              const filter = get(queryFamily(key))
              query[key] = filter
            }
            return query
          },
          set: ({set}, [key, value]) => {
            set(queryFamily(key), isFunction(value)
              ? old => clearEmpty(value(old))
              : clearEmpty(value))
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
    const useResults = () => {
      const setResultsUsed = useSetRecoilState(resultsUsedState)

      // When loading results for the first time, inform the context that
      // results need to be fetched.
      useEffect(() => {
        setResultsUsed(true)
      }, [setResultsUsed])

      return useRecoilValue(resultsState)
    }

    /**
     * Hook for returning an object containing the last used API call.
     *
     * @returns {object} {method, url, body, response}
     */
    const useApiData = () => useRecoilValue(apiDataState)

    /**
     * Hook for modifying an aggregation request and fetching the latest values for
     * this aggregation.
     *
     * @param {string} name The filter name
     * @param {bool} update Whether the hook needs to react to changes in the
     * current query context. E.g. if the component showing the data is not visible,
     * this can be set to false.
     * @param {string} id The aggregation identifier
     * @param {object} config The aggregation configuration. Remember that the
     * value should be a constant or be memoized to avoid a new object being
     * created for each call.
     *
     * @returns {object} An object containing the aggregation results: the layout is
     * specific for each aggregation type.
     */
    const useAgg = (name, update, id = 'default', config = undefined) => {
      const key = useMemo(() => `${name}:${id}`, [name, id])
      const setAgg = useSetRecoilState(aggsFamily(key))
      const aggResponse = useRecoilValue(aggsResponseFamily(key))

      // Whenever the aggregation requirements change, create the final
      // aggregation config and set it in the search context: this will then
      // trigger any required API calls that also return the aggregation
      // response that is returned by this hook.
      useEffect(() => {
        const defaults = filterData[name]?.aggs?.[config?.type]
        const finalConfig = {
          update: update,
          ...defaults || {},
          ...config
        }
        setAgg(finalConfig)
      }, [name, update, setAgg, config])

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
      aggsFamily,
      aggsState,
      paginationState,
      resultsUsedState,
      resultsState,
      apiDataState,
      aggsResponseState,
      isMenuOpenState,
      isCollapsedState,
      isStatisticsEnabledState,
      useStatisticValue,
      useSetStatistic,
      useStatisticState,
      useStatisticsValue,
      useStatisticsState,
      useResults,
      useApiData,
      useAgg,
      useSetFilters
    ]
  }, [initialQuery, filters, filtersLocked, initialStatistics, initialAggs, initialPagination, resource, filterData])

  const setResults = useSetRecoilState(resultsState)
  const setApiData = useSetRecoilState(apiDataState)
  const updateAggsResponse = useSetRecoilState(aggsResponseState)
  const aggs = useRecoilValue(aggsState)
  const query = useRecoilValue(queryState)
  const [pagination, setPagination] = useRecoilState(paginationState)
  const resultsUsed = useRecoilValue(resultsUsedState)
  const updateQueryString = useUpdateQueryString()

  /**
   * This function is used to sync up API calls so that they update the search
   * context state in the same order as they were originally issued.
   *
   * As we cannot guarantee the order in which the API calls finish, we push all
   * calls into a queue. The queue makes sure that API calls get resolved in the
   * original order not matter how long the actual call takes.
   */
  const resolve = useCallback(prop => {
    const {response, timestamp, queryChanged, paginationChanged, search, aggsToUpdate, resource, callback} = prop
    const data = response.response
    const next = apiQueue.current[0]
    if (next !== timestamp) {
      apiMap.current[timestamp] = prop
      return
    }
    // Update the aggregations if new aggregation data is received. The old
    // aggregation data is preserved and new information is updated.
    if (!isEmpty(data.aggregations)) {
      const newAggs = toGUIAgg(data.aggregations, aggsToUpdate, resource)
      callback && callback(newAggs)
      updateAggsResponse(newAggs)
    } else {
      callback && callback(null)
    }
    // Update the query results if new data is received.
    if (queryChanged || paginationChanged) {
      const isExtend = search.pagination.page_after_value
      paginationResponse.current = data.pagination
      setResults(old => {
        const newResults = old ? {...old} : {}
        isExtend ? newResults.data = [...newResults.data, ...data.data] : newResults.data = data.data
        newResults.pagination = combinePagination(search.pagination, data.pagination)
        newResults.setPagination = setPagination
        return newResults
      })
    }
    // Remove this query from queue and see if next can be resolved.
    apiQueue.current.shift()
    setApiData(response)
    const nextTimestamp = apiQueue.current[0]
    const nextResolve = apiMap.current[nextTimestamp]
    if (nextResolve) {
      resolve(nextResolve)
    }
  }, [setApiData, setPagination, setResults, updateAggsResponse])

  /**
   * Function that preprocesses API call requests and finally performs the
   * actual API call.
   *
   * All of the heavier pre-processing, checking, etc. should be done in this
   * function, as it is the final one that gets called after the debounce
   * interval.
   */
  const apiCall = useCallback((query, aggs, pagination, resultsUsed, queryChanged, paginationChanged, updateAggs, callback = undefined) => {
    // Gather aggregations that need to be updated (=require and API call), and
    // the ones that have changed (have changed but do not require an API call).
    const aggsToUpdate = []
    const aggsChanged = []
    for (const [key, agg] of Object.entries(aggs)) {
      if (agg.update) {
        aggsToUpdate.push(key)
      }
      if (agg.changed) {
        aggsChanged.push(key)
      }
    }
    const apiQuery = {...query}
    if (filterDefaults) {
      for (const [key, value] of Object.entries(filterDefaults)) {
        if (isNil(query[key])) {
          apiQuery[key] = value
        }
      }
    }
    const search = {
      owner: apiQuery.visibility,
      query: toAPIFilter(apiQuery, resource),
      aggregations: toAPIAgg(
        aggs,
        resource
      ),
      pagination: {...pagination}
    }

    // When aggregations have changed but the query has not, we request only the
    // aggregation data without any hits. Also if results are not needed, we
    // don't return them.
    if ((updateAggs && !queryChanged && !paginationChanged) || !resultsUsed) {
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

    // When the query changes, information about the most recent agg update (no
    // matter if it causes an API call or not) gets stored. When the query has
    // not changed, the last agg that caused an actual API call is added to the
    // existing list. This way previous aggregation data can be re-used.
    if (queryChanged) {
      updatedAggsMap.current = Object.fromEntries(aggsChanged.map((key) => [key, aggs[key]]))
    } else {
      aggsToUpdate.forEach((key) => { updatedAggsMap.current[key] = aggs[key] })
    }

    // The list of updated filters are always reset.
    updatedFilters.current = new Set()
    firstLoad.current = false
    oldQuery.current = query
    oldPagination.current = pagination

    const timestamp = Date.now()
    apiQueue.current.push(timestamp)
    api.query(resource, search, {loadingIndicator: true, returnRequest: true})
      .then((response) => {
        return resolve({
          response,
          timestamp,
          queryChanged,
          paginationChanged,
          search,
          aggsToUpdate,
          resource,
          callback
        })
      })
      .catch((error) => {
        raiseError(error)
        callback && callback(undefined, error)
      })
  }, [filterDefaults, resource, api, raiseError, resolve])

  // This is a debounced version of apiCall.
  const apiCallDebounced = useMemo(() => debounce(apiCall, 400), [apiCall])

  /**
   * Intermediate function that should primarily be used when trying to perform
   * an API call. Ensures that ensures that:
   * - Calls are debounced when necessary
   * - API calls are made only if necessary
   */
  const apiCallInterMediate = useCallback((query, aggs, pagination, resultsUsed, callback = undefined, forceUpdate = false) => {
    if (disableUpdate.current) {
      disableUpdate.current = false
      return
    }

    // Determine which parts need to be updated. Due to the queueing mechanism
    // we can instantly update the reference to contain the latest information
    // about what will be updated by this query. This ensures that the next query
    // immediately knows the current state even before the API call is finished
    // (or even if no API call is made).
    const queryChanged = forceUpdate ? true : query !== oldQuery.current
    const paginationChanged = pagination !== oldPagination.current
    const [reducedAggs, updateAggs] = reduceAggs(
      aggs,
      updatedAggsMap.current,
      queryChanged,
      updatedFilters.current
    )

    // If results are needed and query and pagination have not changed and
    // aggregations do not need to be updated, no update is necessary. The API
    // calls is made immediately when requesting the first set of results, when
    // the pagination changes or when only aggregations need to be updated.
    // Otherwise it is debounced.
    if ((resultsUsed && (paginationChanged || queryChanged)) || updateAggs) {
      if (firstLoad.current || paginationChanged || !queryChanged) {
        apiCall(query, reducedAggs, pagination, resultsUsed, queryChanged, paginationChanged, updateAggs, callback, false)
      } else {
        apiCallDebounced(query, reducedAggs, pagination, resultsUsed, queryChanged, paginationChanged, updateAggs, callback)
      }
    } else {
      callback && callback(undefined, undefined)
    }
  }, [apiCall, apiCallDebounced])

  // When query, aggregation or pagination changes, update the search context
  useEffect(() => {
    apiCallInterMediate(query, aggs, pagination, resultsUsed)
  }, [query, aggs, pagination, resultsUsed, apiCallInterMediate])

  // This updated the query string to represent the latest value within the
  // search context.
  useEffect(() => {
    updateQueryString()
  }, [updateQueryString])

  // The context contains a set of functions that can be used to hook into
  // different data.
  const values = useMemo(() => {
    // Hook for refreshing the results.
    const useRefresh = () => {
      const query = useRecoilValue(queryState)
      const aggs = useRecoilValue(aggsState)
      const pagination = useRecoilValue(paginationState)
      const resultsUsed = useRecoilValue(resultsUsedState)

      const refresh = useCallback(() => {
        apiCallInterMediate(query, aggs, pagination, resultsUsed, undefined, true)
      }, [aggs, pagination, query, resultsUsed])
      return refresh
    }

    // Hook for imperatively requesting aggregation data. By using this hook you
    // can track the state of individual calls and perform callbacks.
    const useAggCall = (name, id) => {
      const key = useMemo(() => `${name}:${id}`, [name, id])
      const query = useRecoilValue(queryState)
      const pagination = useRecoilValue(paginationState)
      const resultsUsed = useRecoilValue(resultsUsedState)
      const setAgg = useSetRecoilState(aggsFamily(key))

      /**
       * @param {number} size The new aggregation size
       * @param {string} id Identifier for this call
       * @param {function} callback: Function that returns an array containing the
       * new aggregation response and an error if one was encountered. Returns the
       * special value 'undefined' for the response if no update was necessary.
       */
      const aggCall = useCallback((config, callback) => {
        const aggs = {[key]: {...config, update: true}}
        apiCallInterMediate(
          query, aggs, pagination, resultsUsed,
          (response) => callback(response && response[key])
        )

        // We also need to update aggregation request state, otherwise the
        // subsequent calls will not be able to know what was done by this call.
        // To do this without triggering another API call, we disable API updates
        // for one cycle.
        disableUpdate.current = true
        aggs[key].update = false
        setAgg(aggs[key])
      }, [key, query, pagination, resultsUsed, setAgg])

      return aggCall
    }

    return {
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
      useStatisticsState: useStatisticsState,
      useUpdateQueryString: useUpdateQueryString,
      useResults: useResults,
      useApiData: useApiData,
      useAgg: useAgg,
      useAggCall: useAggCall,
      useSetFilters: useSetFilters,
      filters: filters,
      filterData: filterData
    }
  }, [
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
    useStatisticsState,
    useUpdateQueryString,
    useResults,
    useApiData,
    useAgg,
    useSetFilters,
    filters,
    filterData,
    queryState,
    aggsState,
    paginationState,
    resultsUsedState,
    apiCallInterMediate,
    aggsFamily,
    isMenuOpenState,
    isCollapsedState,
    isStatisticsEnabledState
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
  const statistics = {}
  const stats = queryObj.statistics
  if (stats) {
    const statsArray = isArray(stats) ? stats : [stats]
    statsArray.forEach((name, i) => { statistics[name] = {index: i} })
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
      for (const [keyInner, valueInner] of Object.entries(value)) {
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
  // The shown statistics are serialized here: currently we only store a simple
  // ordered list in the query string. If more properties are needed in the
  // future, the whole object should be stored.
  if (!isEmpty(statistics)) {
    queryStringQuery.statistics = Object.entries(statistics)
      .sort((a, b) => a[1].index - b[1].index)
      .map(([key, value]) => key)
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
export function toAPIFilter(query, resource) {
  const queryCustomized = {}
  if (!query) {
    return undefined
  }

  // Perform custom transformations
  function customize(key, value, parent, subKey = undefined) {
    const data = filterDataGlobal[key]

    // Filters that affect the GUI only do not need to be considered
    const guiOnly = data?.guiOnly
    if (guiOnly) {
      return
    }

    // Sections need to be recursively handled. Notice that we cant directly
    // write to an recoil Atom and create a new object for storing the values.
    const section = data?.section
    if (section) {
      const sectionData = {}
      for (const [keyNested, valueNested] of Object.entries(value)) {
        customize(`${key}.${keyNested}`, valueNested, sectionData, keyNested)
      }
      parent[key] = sectionData
    // Regular values are set directly to the parent. A custom setter may be
    // used.
    } else {
      const setter = data?.value?.set
      if (setter) {
        setter(queryCustomized, query, value)
      } else {
        parent[subKey || key] = value
      }
    }
  }
  for (const [k, v] of Object.entries(query)) {
    customize(k, v, queryCustomized)
  }

  // Create the API-compatible keys and values.
  const queryNormalized = {}
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
    const combine = query.combine
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
      newValue = newValue.map((item) => item instanceof Quantity
        ? item.toSI().value()
        : item)
    }
  } else if (value instanceof Quantity) {
    newValue = value.toSI().value()
  } else if (isArray(value)) {
    if (value.length === 0) {
      newValue = undefined
    } else {
      newValue = value.map((item) => item instanceof Quantity
        ? item.toSI().value()
        : item)
    }
  } else if (isPlainObject(value)) {
    newValue = {}
    for (const [keyInner, valueInner] of Object.entries(value)) {
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
    for (const [key, value] of Object.entries(query)) {
      const newKey = filterFullnames[key] || key
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
    for (const [keyInner, valueInner] of Object.entries(value)) {
      const valueConverted = toGUIFilterSingle(keyInner, valueInner, units, fullPath)
      if (!isNil(valueConverted)) {
        newValue[keyInner] = valueConverted
      }
    }
  } else {
    // If the key is an operator, the filter name is read from the path.
    const opKeys = new Set(['lte', 'lt', 'gte', 'gt'])
    const filterPath = opKeys.has(key) ? path : fullPath
    const multiple = filterDataGlobal[filterPath].multiple
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
 * @param {object} updatedFilters Set of filters that were updated together with
 * this call.
 * @param {string} resource The resource we are looking at: entries or materials.
 *
 * @returns {object} Aggregation query that is usable by the API.
 */
function toAPIAgg(aggs, resource) {
  const apiAggs = {}
  for (const [key, agg] of Object.entries(aggs)) {
    const filter_name = key.split(':')[0]
    if (agg.update) {
      const exclusive = filterDataGlobal[filter_name].exclusive
      const type = agg.type
      const name = resource === 'materials' ? materialNames[filter_name] : filter_name
      const apiAgg = apiAggs[name] || {}
      const aggSet = agg.set
      const finalAgg = aggSet ? aggSet(agg) : agg
      apiAgg[type] = {
        quantity: name,
        // Exclusive quantities (quantities that have one value per entry) are
        // always fetched with exclude_from_search
        exclude_from_search: exclusive,
        ...finalAgg
      }
      apiAggs[key] = apiAgg
    }
  }
  return apiAggs
}

/**
 * Used to transform an API aggregation result into a form that is usable by the
 * GUI.
 *
 * @param {object} aggs The aggregation data as returned by the API.
 * @param {array} aggsToUpdate The set of targeted filters. Needed because the keys
 * in the aggs dictionary may be different due to custom aggregation set/get.
 * @param {string} resource The resource we are looking at: entries or materials.
 *
 * @returns {object} Aggregation result that is usable by the GUI.
 */
function toGUIAgg(aggs, aggsToUpdate, resource) {
  if (isEmpty(aggs)) {
    return {}
  }
  // Modify keys according to target resource (entries/materials).
  let aggsNormalized
  if (resource === 'materials') {
    aggsNormalized = {}
    for (const key of Object.keys(aggs)) {
      const filter_name = key.split(':')[0]
      const name = resource === 'materials' ? entryNames[filter_name] : filter_name
      aggs[key].quantity = name
      aggsNormalized[key] = aggs[key]
    }
  } else {
    aggsNormalized = aggs
  }

  // Perform custom transformations
  const aggsCustomized = {}
  for (const key of aggsToUpdate) {
    const filter_name = key.split(':')[0]
    aggs = aggsNormalized[key]
    if (!isNil(aggs)) {
      for (const [type, agg] of Object.entries(aggs)) {
        const aggGet = filterDataGlobal[filter_name]?.aggs?.[type].get
        const aggFinal = aggGet ? aggGet(agg) : agg
        // Add flag for if all terms have been returned, and the total number of
        // items. TODO: Could this total be given by the API directly?
        if (type === 'terms') {
          aggFinal.exhausted = aggFinal.size > aggFinal.data.length
          aggFinal.total = aggFinal.data.reduce((a, b) => a + b.count, 0)
        }
        aggsCustomized[key] = aggFinal
      }
    }
  }
  return aggsCustomized
}

/**
 * Goes through the given aggregations and compares then against the old results
 * and the current query to identify and return the aggregations that need to be
 * called.
 *
 * @param {object} aggs The current aggregation configuration.
 * @param {object} oldAggs Reduced aggregation config from latest finished query.
 * @param {bool} queryChanged Whether the query has changed.
 *
 * @returns {object} Reduced aggregation config.
 */
function reduceAggs(aggs, oldAggs, queryChanged, updatedFilters) {
  // Loop through the different aggregations for each quantity and see if any
  // of them need to be updated.
  const reducedAggs = {}
  let updateAggs = false
  for (const [key, agg] of Object.entries(aggs)) {
    const filter_name = key.split(':')[0]
    if (!isBoolean(agg.update)) {
      throw Error(`It was not specified whether the aggregation ${key} should update or not.`)
    }
    let update = agg.update
    let changed = agg.update

    // If the old aggregation data is a superset of the newly queried
    // aggregation data, and the query has not changed, there is no need to
    // update.
    const oldAgg = oldAggs?.[key]
    if (update && !queryChanged && oldAgg) {
      const aggParams = {...agg}
      const aggParamsOld = {...oldAgg}
      delete aggParams.changed
      delete aggParams.update
      delete aggParamsOld.changed
      delete aggParamsOld.update
      if (isEqual(aggParams, aggParamsOld)) {
        update = false
        changed = false
      } else if (!isNil(agg.size)) {
        delete aggParams.size
        delete aggParamsOld.size
        if (isEqual(aggParams, aggParamsOld) && oldAgg.size >= agg.size) {
          update = false
        }
      }
    }

    // If the filter is exclusive, and ONLY it has been modified in this
    // query, we do not update it's aggregation.
    const exclude = isNil(agg.exclude_from_search)
      ? filterDataGlobal[filter_name].exclusive
      : agg.exclude_from_search
    if (exclude && updatedFilters.has(filter_name) && updatedFilters.size === 1) {
      update = false
    }

    // Create the new reduced version
    const newAgg = {...agg, update, changed}
    if (update) {
      updateAggs = true
    }
    reducedAggs[key] = newAgg
  }
  return [reducedAggs, updateAggs]
}
