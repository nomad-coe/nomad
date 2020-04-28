import React, { useState, useContext, useEffect, useRef, useCallback } from 'react'
import PropTypes from 'prop-types'
import hash from 'object-hash'
import { errorContext } from '../errors'
import { onlyUnique, objectFilter } from '../../utils'
import { domains } from '../domains'
import { apiContext } from '../api'
import { useLocation, useHistory } from 'react-router-dom'
import qs from 'qs'
import * as searchQuantities from '../../searchQuantities.json'

const padDateNumber = number => String('00' + number).slice(-2)

export const Dates = {
  dateHistogramStartDate: '2014-12-15',
  APIDate: date => date.toISOString(),
  JSDate: date => new Date(date),
  FormDate: date => {
    date = new Date(date)
    return `${date.getFullYear()}-${padDateNumber(date.getMonth())}-${padDateNumber(date.getDate())}`
  },
  addSeconds: (date, interval) => new Date(date.getTime() + interval * 1000),
  deltaSeconds: (from, end) => Math.round((new Date(end).getTime() - new Date(from).getTime()) / 1000),
  intervalSeconds: (from, end, buckets) => Math.round((new Date(end).getTime() - new Date(from).getTime()) / (1000 * buckets)),
  buckets: 50
}

searchQuantities['from_time'] = true
searchQuantities['until_time'] = true
/**
 * A custom hook that reads and writes search parameters from the current URL.
 */
const useSearchUrlQuery = () => {
  const location = useLocation()
  const history = useHistory()
  const urlQuery = location.search ? {
    ...qs.parse(location.search.substring(1))
  } : {}
  const searchQuery = objectFilter(urlQuery, key => searchQuantities[key])
  const rest = objectFilter(urlQuery, key => !searchQuantities[key])
  if (searchQuery.atoms && !Array.isArray(searchQuery.atoms)) {
    searchQuery.atoms = [searchQuery.atoms]
  }
  if (searchQuery.only_atoms && !Array.isArray(searchQuery.only_atoms)) {
    searchQuery.only_atoms = [searchQuery.only_atoms]
  }
  const setUrlQuery = query => {
    history.push(location.pathname + '?' + qs.stringify({
      ...rest,
      ...objectFilter(query, key => query[key])
    }, {indices: false}))
  }
  return [searchQuery, setUrlQuery]
}

/**
 * The React context object. Can be accessed from functional components with useContext.
 */
export const searchContext = React.createContext()

/**
 * Component that provides a searchContext. Can be used with useContext. The context
 * objects provides access to the current search request and response as well as
 * callbacks to manipulate the current search request.
 *
 * The search request is made from two objects: the request and the query. The former
 * contains all parameters that do not effect the search results themselves. This includes
 * pagination, statistics, order. The query object contains all parameters that
 * constitute the actual search. This includes the domain and owner parameters.
 */
export default function SearchContext({initialRequest, initialQuery, query, children}) {
  const defaultStatistics = [] // ['atoms', 'authors'] TODO
  const emptyResponse = {
    statistics: {
      total: {
        all: {}
      }
    },
    pagination: {
      total: undefined,
      per_page: 10,
      page: 1,
      order: -1,
      order_by: 'upload_time'
    },
    metric: domains.dft.defaultSearchMetric
  }

  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)

  // React calls the children effects for the parent effect. But the parent effect is
  // run with the state of the last render, which is the state before the children effects.
  // If we would maintain the request in regular React state, we might execute unnecessary
  // outdated requests.
  // Therefore, we use two ref objects and one state object to manage the current state.
  // The goal is to reduce the amounts of re-renders, only send requests to the api with
  // the latest set parameters, and only send requests if necessary.
  // The first ref keeps all information that will form the next
  // search request that is send to the API. It therefore keeps the state of the current
  // request provided by the various children of this context. It also helps us to
  // lower the amount of state changes.
  // The second ref keeps a hash over the last request that was send to the API.
  // This is used to verify if a new request is actually necessary.
  const requestRef = useRef({
    metric: domains.dft.defaultSearchMetric,
    statistics: [],
    groups: {},
    domainKey: domains.dft.key,
    owner: 'all',
    pagination: {
      page: 1,
      per_page: 10,
      order: -1,
      order_by: 'upload_time'
    },
    statisticsToRefresh: [],
    query: {},
    update: 0
  })
  const lastRequestHashRef = useRef(0)

  // We use proper React state to maintain the last response from the API.
  const [response, setResponse] = useState(emptyResponse)

  // We use a custom hook to read/write search parameters from the current URL.
  const [urlQuery, setUrlQuery] = useSearchUrlQuery()
  // We set the current query. This will be used by an effect to potentially call the
  // API after rendering.
  requestRef.current.query = urlQuery

  // This is a callback that executes the current request in requestRef without any
  // checks for necessity. It will update the response state, once the request has
  // been answered by the api.
  const runRequest = useCallback(() => {
    let dateHistogramInterval = null
    const {metric, domainKey, owner, dateHistogram} = requestRef.current
    const domain = domains[domainKey]
    const apiRequest = {
      ...initialRequest,
      ...requestRef.current.pagination,
      statistics: requestRef.current.statistics,
      ...requestRef.current.groups,
      metrics: (metric === domain.defaultSearchMetric) ? [] : [metric],
      domain: domain.key
    }
    const apiQuery = {
      ...apiRequest,
      owner: owner,
      ...initialQuery,
      ...requestRef.current.query,
      ...query
    }
    if (dateHistogram) {
      dateHistogramInterval = Dates.intervalSeconds(
        apiQuery.from_time || Dates.dateHistogramStartDate,
        apiQuery.until_time || new Date(), Dates.buckets)
      apiQuery['date_histogram'] = true
      apiQuery['interval'] = `${dateHistogramInterval}s`
    }
    api.search(apiQuery)
      .then(newResponse => {
        setResponse({
          ...emptyResponse,
          ...newResponse,
          metric: metric,
          dateHistogramInterval: dateHistogramInterval,
          from_time: apiQuery.from_time,
          until_time: apiQuery.until_time
        })
      }).catch(error => {
        setResponse({...emptyResponse, metric: metric})
        raiseError(error)
      })
  }, [requestRef, setResponse, api])

  // The following are various callbacks that can be used by children to update the
  // request and implicitly trigger a search request to the API. The implicit triggering
  // is realised that all changes to the request are accompanied by updates to the URL
  // which is used to hold the whole request state. Each push to the history will rerender
  // everything and therefore trigger effects.
  const setRequestParameters = useCallback(
    changes => {
      requestRef.current.pagination = {
        ...requestRef.current.pagination,
        ...changes
      }
    }, [requestRef])

  const setDomain = useCallback(domainKey => {
    requestRef.current.domainKey = domainKey || domains.dft.key
  }, [requestRef])

  const setOwner = useCallback(owner => {
    requestRef.current.owner = owner
  }, [requestRef])

  const setMetric = useCallback(metric => {
    requestRef.current.metric = metric || domains.dft.defaultSearchMetric
  }, [requestRef])

  const setStatistics = useCallback(statistics => {
    requestRef.current.statistics = [...statistics, ...defaultStatistics].filter(onlyUnique)
  }, [requestRef])

  const setGroups = useCallback(groups => {
    requestRef.current.groups = {...groups}
  }, [requestRef])

  const setDateHistogram = useCallback(dateHistogram => {
    requestRef.current.dateHistogram = dateHistogram
  }, [requestRef])

  const handleQueryChange = (changes, replace) => {
    if (changes.atoms && changes.atoms.length === 0) {
      changes.atoms = undefined
    }
    if (changes.only_atoms && changes.only_atoms.length === 0) {
      changes.only_atoms = undefined
    }

    if (replace) {
      setUrlQuery(changes)
    } else {
      setUrlQuery({...urlQuery, ...changes})
    }
  }

  // We check and run (if necessary) the search request after each render
  useEffect(() => {
    if (lastRequestHashRef.current !== hash(requestRef.current)) {
      runRequest()
      lastRequestHashRef.current = hash(requestRef.current)
    }
  })

  const value = {
    response: response,
    query: {
      ...requestRef.current.query
    },
    apiQuery: {
      domain: requestRef.current.domainKey,
      owner: requestRef.current.owner,
      ...requestRef.current.query,
      ...query
    },
    domain: domains[requestRef.current.domainKey],
    metric: requestRef.current.metric,
    requestParameters: requestRef.current.pagination,
    setRequestParameters: setRequestParameters,
    setQuery: handleQueryChange,
    setMetric: setMetric,
    setGroups: setGroups,
    setDomain: setDomain,
    setOwner: setOwner,
    setStatisticsToRefresh: () => null, // TODO remove
    setStatistics: setStatistics,
    setDateHistogram: setDateHistogram,
    update: runRequest
  }

  return <searchContext.Provider value={value} >{children}</searchContext.Provider>
}
SearchContext.propTypes = {
  /**
   * An object with initial query parameters. These will be added to the search context
   * and be used in all search requests.
  */
  query: PropTypes.object,
  /**
   * An object with initial query parameters. These will be added to the search context
   * and the first search request. Afterwards search parameters might be removed or
   * overwritten by the search.
   */
  initialQuery: PropTypes.object,
  /**
   * An object with initial request parameters. These will be added to the search context
   * and the first search request. Afterwards request parameters might be removed or
   * overwritten by children components.
   */
  initialRequest: PropTypes.object,
  /**
   * The children prop. All components in the children can make use of this search
   * context via useContext.
   */
  children: PropTypes.any
}
