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
import { atom, useSetRecoilState, useRecoilValue } from 'recoil'
import _ from 'lodash'
// import { useHistory } from 'react-router-dom'
// import qs from 'qs'
import { useApi } from '../apiV1'
import { setToArray } from '../../utils'

// This holds the shared query that is controlled by adding/modifying/removing
// filters.
export const queryStringState = atom({
  key: 'queryString',
  default: {}
})
export const queryState = atom({
  key: 'query',
  default: {}
})
export const loadingState = atom({
  key: 'loading',
  default: false
})

export function useLoadingValue() {
  const loading = useRecoilValue(loadingState)
  return loading
}

// This holds the shared statistics that is controllled by
// adding/modifying/removing filters.
export const statisticsState = atom({
  key: 'statistics',
  default: []
})

// This holds the shared statistics that is controllled by
// adding/modifying/removing filters.
export const aggregationsState = atom({
  key: 'aggregations',
  default: []
})

export function useQueryState() {
  const setQuery = useSetQuery()
  const query = useQueryValue()

  return [query, setQuery]
}

export function useSetQuery() {
  const setQuery = useSetRecoilState(queryState)
  // const setQueryString = useSetRecoilState(queryStringState)
  // const history = useHistory()

  // function setQuery(func) {
  //   setQueryString(oldQuery => {
  //     const newQuery = qs.parse(oldQuery)
  //     const temp = func(newQuery)
  //     const elements = temp['results.material.elements']
  //     if (elements) {
  //       temp['results.material.elements'] = [...elements]
  //     }
  //     console.log(temp)
  //     const queryString = qs.stringify(temp, {indices: false})
  //     // history.push(history.location.pathname + '?' + queryString)
  //     return queryString
  //   })
  // }
  return setQuery
}

export function useQueryValue() {
  const query = useRecoilValue(queryState)
  // const queryString = useRecoilValue(queryStringState)
  // const temp = qs.parse(queryString)
  // let elements = temp['results.material.elements']
  // if (elements) {
  //   if (!Array.isArray(elements)) {
  //     elements = [elements]
  //   }
  //   temp['results.material.elements'] = new Set(elements)
  // }
  return query
}

/**
 * This hook will expose a function for setting filter values for a specific
 * quantity. Use this hook if you intend to only set the filter value and are
 * not interested in the query results.
 *
 * @param {*} quantity Name of the quantity to set. Should exist in searchQuantities.json.
 * @returns function for setting the value for the given quantity
 */
export function useSetFilter(quantity, set) {
  set.add(quantity)
  const setQuery = useSetQuery()

  function setFilter(value) {
    setQuery(oldQuery => {
      const newQuery = {...oldQuery}
      newQuery[quantity] = value
      return newQuery
    })
  }

  return setFilter
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
  const query = useQueryValue()
  return query[quantity]
}

/**
 * This hook will expose a function for getting and setting filter values for a
 * specific quantity. Use this hook if you intend to both read and write the
 * filter value.
 *
 * @param {*} quantity Name of the quantity to set. Should exist in searchQuantities.json.
 * @returns array containing the filter value and setter function for it.
 */
export function useFilterState(quantity, set) {
  const setFilter = useSetFilter(quantity, set)
  const filter = useFilterValue(quantity)

  return {
    filter: filter,
    setFilter: setFilter
  }
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
  const [results, setResults] = useState()
  const setLoading = useSetRecoilState(loadingState)

  // The results are fetched as a side effect in order to not block the
  // rendering. This causes two renders: first one without the data, the second
  // one with the data.
  const apiCall = useCallback(search => {
    const finalSearch = {...search}
    const finalQuery = {...search.query}
    finalQuery['results.material.elements'] = setToArray(finalQuery['results.material.elements'])
    finalSearch.query = finalQuery
    api.queryEntry(finalSearch)
      .then(data => setResults(data))
      .finally(() => setLoading(false))
  }, [api])

  // This is a debounced version of apiCall.
  const debounced = useCallback(_.debounce(apiCall, delay), [])

  // The API call is made immediately on first render. On subsequent renders it
  // will be debounced.
  useEffect(() => {
    if (firstRender.current) {
      apiCall(search)
      firstRender.current = false
    } else {
      setLoading(true)
      debounced(search)
    }
  }, [apiCall, debounced, search])

  return {
    results: results,
    search: search
  }
}

/**
 */
export function useSearch() {
  const stats = useRecoilValue(statisticsState)
  const aggs = useRecoilValue(aggregationsState)
  const query = useQueryValue()
  const result = useMemo(() => {
    return {
      owner: 'all',
      query: query,
      pagination: {
        page: 1,
        page_size: 10,
        order: 'desc',
        order_by: 'upload_time'
      },
      stats: stats,
      aggregations: aggs
    }
  }, [query, stats, aggs])
  return result
}
