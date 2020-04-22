import React, { useState, useContext, useEffect, useRef, useCallback } from 'react'
import PropTypes from 'prop-types'
import hash from 'object-hash'
import { errorContext } from '../errors'
import { onlyUnique } from '../../utils'
import { domains } from '../domains'
import { apiContext } from '../api'

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
export default function SearchContext({initialRequest, initialQuery, children}) {
  const defaultStatistics = ['atoms', 'authors']
  const emptyResponse = {
    statistics: {
      total: {
        all: {}
      }
    },
    pagination: {
      total: undefined,
      per_page: 10,
      page: 1
    },
    metric: domains.dft.defaultSearchMetric
  }

  const {api, info} = useContext(apiContext)
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
  // The state requestHash is used to trigger the effect that will execute the request
  // if necessary. Thus any requests are only send by effects on this context component.
  // If we would send the requests in child effects, we would send unnecessary requests
  // if two children change the request in the same render cycle.
  const requestRef = useRef({
    metric: domains.dft.defaultSearchMetric,
    statistics: [],
    groups: {},
    domainKey: domains.dft.key,
    pagination: {
      order_by: 'upload_id',
      order: -1,
      page: 1,
      per_page: 10
    },
    query: {
      owner: 'all'
    },
    update: 0
  })
  const requestHashRef = useRef(0)
  const [requestHash, setRequestHash] = useState(0)

  // We use proper React state to maintain the last response from the API.
  const [response, setResponse] = useState(emptyResponse)
  const [statisticsToRefresh, setStatisticsToRefresh] = useState([]) // TODO

  // This is a callback that executes the current request in requestRef without any
  // checks for necessity. It will update the response state, once the request has
  // been answered by the api.
  const runRequest = useCallback(() => {
    const {metric, domainKey} = requestRef.current
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
      ...initialQuery,
      ...requestRef.current.query
    }
    api.search(apiQuery, statisticsToRefresh)
      .then(newResponse => {
        setResponse({...emptyResponse, ...newResponse, metric: metric})
      }).catch(error => {
        setResponse({...emptyResponse, metric: metric})
        raiseError(error)
      })
  }, [requestRef, setResponse, statisticsToRefresh, api])

  // This callback will update the requestHash state. This will trigger an effect if
  // the new hash is different from the last.
  // This callback should be called after the requestRef was changed.
  const onRequestChange = useCallback(
    () => {
      setRequestHash(hash(requestRef.current))
    }, [requestRef, setRequestHash]
  )

  // This callback increased the update counter in requestRef therefore causes to re-run
  // the request, even if no parameter have changed. This can be used by children to
  // refresh the search results.
  const update = useCallback(() => {
    requestRef.current.update = requestRef.current.update + 1
    onRequestChange()
  }, [onRequestChange, requestRef])

  // The following are various callbacks that can be used by children to update the
  // request and implicitly trigger a search request to the API.
  const setRequestParameters = useCallback(
    changes => {
      requestRef.current.pagination = {
        ...requestRef.current.pagination,
        ...changes
      }
      onRequestChange()
    }, [onRequestChange, requestRef])

  const setDomain = useCallback(domainKey => {
    requestRef.current.domainKey = domainKey || domains.dft.key
    onRequestChange()
  }, [onRequestChange, requestRef])

  const setOwner = useCallback(owner => {
    requestRef.current.query.owner = owner
    onRequestChange()
  }, [onRequestChange, requestRef])

  const setMetric = useCallback(metric => {
    requestRef.current.metric = metric || domains.dft.defaultSearchMetric
    onRequestChange()
  }, [onRequestChange, requestRef])

  const setStatistics = useCallback(statistics => {
    requestRef.current.statistics = [...statistics, ...defaultStatistics].filter(onlyUnique)
    onRequestChange()
  }, [onRequestChange, requestRef])

  const setGroups = useCallback(groups => {
    requestRef.current.groups = {...groups}
    onRequestChange()
  }, [onRequestChange, requestRef])

  const handleStatisticsToRefreshChange = statistics => setStatisticsToRefresh(
    [...statisticsToRefresh, statistics].filter(onlyUnique)
  )
  const handleQueryChange = (changes, replace) => {
    if (changes.atoms && changes.atoms.length === 0) {
      changes.atoms = undefined
    }
    if (changes.only_atoms && changes.only_atoms.length === 0) {
      changes.only_atoms = undefined
    }

    if (replace) {
      requestRef.current.query = {...changes}
    } else {
      requestRef.current.query = {...requestRef.current.query, ...changes}
    }
    onRequestChange()
  }

  // We initially trigger a search request on mount.
  useEffect(() => {
    // In some cases, especially on mount, requestHash might not be based on the
    // most current requestRef and we have to recompute the hash
    const newRequestHash = hash(requestRef.current)
    if (requestHashRef.current !== newRequestHash) {
      requestHashRef.current = newRequestHash
      runRequest()
    }
  }, [requestHashRef, requestHash, runRequest])

  const value = {
    response: response,
    query: {
      domain: requestRef.current.domainKey,
      ...requestRef.current.query
    },
    domain: domains[requestRef.current.domainKey],
    metric: requestRef.current.metric,
    setRequestParameters: setRequestParameters,
    setQuery: handleQueryChange,
    setMetric: setMetric,
    setGroups: setGroups,
    setDomain: setDomain,
    setOwner: setOwner,
    setStatisticsToRefresh: handleStatisticsToRefreshChange,
    setStatistics: setStatistics,
    update: update
  }

  return <searchContext.Provider value={value} >{children}</searchContext.Provider>
}
SearchContext.propTypes = {
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
