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
import compose from 'compose-function'
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
  isBoolean,
  isEmpty,
  isEqual,
  isArray,
  isPlainObject,
  isNil,
  isSet,
  isFunction,
  size,
  cloneDeep
} from 'lodash'
import qs from 'qs'
import { useLocation, useHistory } from 'react-router-dom'
import { Box, Link, makeStyles } from '@material-ui/core'
import { v4 as uuidv4 } from 'uuid'
import PropTypes from 'prop-types'
import { useApi } from '../api'
import {
  authorList,
  entryName,
  entryType,
  formatTimestamp,
  getDeep,
  formatNumber,
  setToArray,
  parseQuantityName,
  rsplit,
  parseOperator,
  getSuggestions
} from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { useErrors } from '../errors'
import { combinePagination, addColumnDefaults } from '../datatable/Datatable'
import UploadStatusIcon from '../uploads/UploadStatusIcon'
import { getWidgetsObject } from './widgets/Widget'
import { inputSectionContext } from './input/InputSection'
import { SearchSuggestion } from './SearchSuggestion'
import { withFilters } from './FilterRegistry'
import { useUnitContext } from '../units/UnitContext'

const useWidthConstrainedStyles = makeStyles(theme => ({
  root: {
    maxWidth: '500px',
    textOverflow: 'ellipsis',
    overflow: 'hidden'
  }
}))

function WidthConstrained(props) {
  const styles = useWidthConstrainedStyles()
  return <Box className={styles.root}>
    {props.children}
  </Box>
}
WidthConstrained.propTypes = {
  children: PropTypes.node
}

const debounceTime = 450

function parseQueryString(queryString) {
  return queryString
    ? qs.parse(queryString, {comma: true})
    : undefined
}

/**
 * Removes keycloak-specific fragments from the given URL.
*/
function removeKeycloak(url) {
  for (const pattern of ['#state=', '#iss=', '#error=']) {
    url = url.split(pattern)[0]
  }
  return url
}

const RouteListener = React.memo(() => {
  const location = useLocation()
  const history = useHistory()
  const { useQueryString, useSetFilters, useParseQueries } = useSearchContext()
  const queryString = useQueryString()()
  const setFilters = useSetFilters()
  const parseQueries = useParseQueries()
  const oldQuery = useRef()

  // When query changes, update history
  useEffect(() => {
    history.replace(history.location.pathname + '?' + queryString)
    oldQuery.current = queryString
  }, [history, queryString])

  // When location changes, update SearchContext if query has changed.
  useEffect(() => {
    // The query string needs to contain everything in the hash as well, since
    // the syntax for data within custom schemas uses a hashtag. TODO: Keycloak
    // seems to use the hashtag for special purposes as well, which causes some
    // clashes when loading the page from scratch.
    const queryString = removeKeycloak(location.search.slice(1) + location.hash)
    // Empty query string preserves the old filters, as this will enable people
    // to preserve filters when changing apps. This causes minor
    // inconsistencies, but might be better for the UX.
    if (!isEmpty(queryString) && queryString !== oldQuery.current) {
      setFilters(parseQueries(parseQueryString(queryString)))
    }
  }, [location, setFilters, parseQueries])

  return null
})

RouteListener.propTypes = {
}

/**
 * HOC that handles the two-way syncing of URL query parameters and
 * SearchContext filter values.
 */
export const withQueryString = (WrappedComponent) => {
  const WithQueryString = ({initialFilterValues, children, ...rest}) => {
    const firstRender = useRef(true)
    const valueRef = useRef()

    // Extract URL query parameters, ignore keycloak authentication state
    // parameters starting with '#iss='.
    useMemo(() => {
      if (firstRender.current) {
        const split = window.location.href.split('?')
        const search = (split.length > 1)
          ? parseQueryString(removeKeycloak(split.slice(-1)[0]))
          : {}
        valueRef.current = {...(initialFilterValues || {}), ...search}
        firstRender.current = false
      }
    }, [initialFilterValues])

    return <WrappedComponent {...rest} initialFilterValues={valueRef.current}>
      {children}
      <RouteListener />
    </WrappedComponent>
  }

  WithQueryString.displayName = `withFilter(${WrappedComponent.displayName || WrappedComponent.name})`
  WithQueryString.propTypes = {
    initialFilterValues: PropTypes.object, // Determines which filters are available
    children: PropTypes.node // Determines which filters are available
  }

  return WithQueryString
}

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
export const searchContext = React.createContext()
export const SearchContextRaw = React.memo(({
  resource,
  initialFiltersLocked,
  initialColumns,
  initialRows,
  initialFilterMenus,
  initialPagination,
  initialDashboard,
  initialFilterData,
  initialFilterGroups,
  initialFilterValues,
  initialSearchSyntaxes,
  id,
  suggestionHistorySize,
  children
}) => {
  const {api, user} = useApi()
  const {units} = useUnitContext()
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
  const indexFilters = useRef(0)
  const indexLocked = useRef(0)
  const contextID = useMemo(() => uuidv4(), [])

  // Initialize the set of filters that are available in this context
  const {initialFilterPaths, initialFilterAbbreviations} = useMemo(() => {
    return {
      initialFilterPaths: new Set(Object.keys(initialFilterData)),
      initialFilterAbbreviations: getAbbreviations(initialFilterData)
    }
  }, [initialFilterData])

  // The final set of columns
  const columns = useMemo(() => {
    if (!initialColumns) return undefined

    // Add unit information if one is defined. This unit is currently fixed and
    // not affected by global unit system.
    const config = cloneDeep(initialColumns)
    let options = config?.options ? Object.values(config.options) : []
    options.forEach(option => {
      const unit = option.unit
      option.label = option.label || initialFilterData[option.quantity || option.key]?.label
      if (unit) {
        option.unit = new Unit(unit)
        option.label = `${option.label} (${option.unit.label()})`
      }
    })

    // Automatically determine the render function based on metainfo.
    options.forEach(option => {
      option.render = (data) => {
        const value = getDeep(data, parseQuantityName(option.key).path)
        if (isNil(value)) return value

        const transform = (value) => {
          const key = option.key
          const unit = option.unit
          const format = option.format
          if (unit) {
            const originalUnit = initialFilterData[key]?.unit
            value = new Quantity(value, originalUnit).to(unit).value()
          }
          if (format) {
            const dtype = initialFilterData[key]?.dtype
            value = formatNumber(value, dtype, format?.mode, format?.decimals)
          }
          return value
        }
        return isArray(value)
          ? value.map(transform).join(option?.format?.separator || ', ')
          : transform(value)
      }
    })

    // Custom render and other setting overrides are used for a subset of
    // columns.
    const overrides = {
      entry_name: {
        render: (entry) => <WidthConstrained>{entryName(entry)}</WidthConstrained>
      },
      entry_type: {
        render: (entry) => <WidthConstrained>{entryType(entry)}</WidthConstrained>
      },
      upload_create_time: {
        render: row => row?.upload_create_time
          ? formatTimestamp(row.upload_create_time)
          : <i>no upload time</i>
      },
      entry_create_time: {
        render: row => row?.entry_create_time
          ? formatTimestamp(row.entry_create_time)
          : <i>no entry time</i>
      },
      authors: {
        sortable: false,
        render: row => authorList(row),
        align: 'left'
      },
      references: {
        sortable: false,
        render: row => {
          const refs = row.references || []
          if (refs.length > 0) {
            return (
              <div style={{display: 'inline'}}>
                {refs.map((ref, i) => <span key={ref}>
                  <Link href={ref}>{ref}</Link>{(i + 1) < refs.length ? ', ' : <React.Fragment/>}
                </span>)}
              </div>
            )
          } else {
            return <i>no references</i>
          }
        }
      },
      datasets: {
        sortable: false,
        render: entry => {
          const datasets = entry.datasets || []
          if (datasets.length > 0) {
            return datasets.map(dataset => dataset.dataset_name).join(', ')
          } else {
            return <i>no datasets</i>
          }
        }
      },
      published: {
        render: (entry) => <UploadStatusIcon data={entry} user={user} />
      }
    }

    // Sort options by putting initially selected ones on top
    if (config.selected) {
      options = [
        ...config.selected.map(key => options.find(opt => opt.key === key)).filter(opt => !!opt),
        ...options.filter(opt => !config.selected.find(key => key === opt.key))
      ]
    }

    // Add defaults and custom overrides to the options
    addColumnDefaults(options, undefined, initialFilterData)
    config.options = Object.fromEntries(options.map(option => {
      return [option.key, {...option, ...(overrides[option.key] || {})}]
    }))

    return config
  }, [initialFilterData, initialColumns, user])

  // The final row configuration
  const rows = useMemo(() => {
    return initialRows || undefined
  }, [initialRows])

  // The final menu configuration
  const filterMenus = useMemo(() => {
    return initialFilterMenus || undefined
  }, [initialFilterMenus])

  // The final dashboard configuration
  const dashboard = useMemo(() => {
    return initialDashboard
      ? {...initialDashboard, widgets: getWidgetsObject(initialDashboard?.widgets || [])}
      : undefined
  }, [initialDashboard])

  // Initialize the search context state: any parameters in the URL are read and
  // default values as specified in filter registry are loaded
  const [initialQuery, initialAggs, filterDefaults] = useMemo(() => {
    const filterDefaults = {}
    for (const [key, value] of Object.entries(initialFilterData)) {
      if (!isNil(value.default)) {
        filterDefaults[key] = value.default
      }
    }
    const initialAggs = {}
    for (const key of initialFilterPaths) {
      initialAggs[key] = {
        default: {update: false}
      }
    }
    return [
      parseQueries(initialFilterValues, initialFilterData, initialFilterAbbreviations.filterFullnames),
      initialAggs,
      filterDefaults
    ]
  }, [initialFilterData, initialFilterAbbreviations, initialFilterValues, initialFilterPaths])

  // Atoms + setters and getters are used instead of regular React states to
  // avoid re-rendering components that are not depending on these values. The
  // Recoil mechanisms are abstracted under hooks in order to make it easier to
  // replace Recoil with some other state management library in the future.
 const {
      useFilterLocked,
      useFiltersLocked,
      useFiltersLockedState,
      useFilterValue,
      useSetFilter,
      useFilterState,
      useFilters,
      useFilterNames,
      useFiltersData,
      useFilterMaps,
      useResetFilters,
      useQueryString,
      useDynamicQueryModes,
      aggsFamily,
      useWidgetValue,
      useSetWidget,
      useWidgetState,
      useWidgetsValue,
      useWidgetsState,
      useAddWidget,
      useRemoveWidget,
      useResetWidgets,
      useQuery,
      useSearchSuggestionsValue,
      useSetSearchSuggestions,
      usePagination,
      useSetPagination,
      usePaginationState,
      useAgg,
      useAggs,
      useSetAggsResponse,
      useSetFilters,
      useUpdateFilter,
      useIsMenuOpen,
      useSetIsMenuOpen,
      useIsCollapsed,
      useSetIsCollapsed,
      useIsStatisticsEnabled,
      useSetIsStatisticsEnabled,
      useApiQuery,
      useSetApiQuery,
      useApiData,
      useSetApiData,
      useResults,
      useSetResults,
      useResultsUsed,
      useRequired
  } = useMemo(() => {
    const isMenuOpenState = atom({
      key: `isMenuOpen_${contextID}`,
      default: false
    })

    const isCollapsedState = atom({
      key: `isCollapsed_${contextID}`,
      default: false
    })

    const isStatisticsEnabledState = atom({
      key: `statisticsEnabled_${contextID}`,
      default: true
    })

    const apiDataState = atom({
      key: `apiData_${contextID}`,
      default: null
    })

    const apiQueryState = atom({
      key: `apiQuery_${contextID}`,
      default: null
    })

    const resultsState = atom({
      key: `results_${contextID}`,
      default: {
        pagination: {}
      }
    })

    const resultsUsedState = atom({
      key: `resultsUsed_${contextID}`,
      default: false
    })

    let suggestions = []
    const suggestionID = id && `nomad-searchcontext-${id}`
    if (suggestionID) {
      try {
        suggestions = JSON
          .parse(window.localStorage.getItem(suggestionID))
          .map(item => new SearchSuggestion(item))
          .slice(0, suggestionHistorySize)
      } catch (e) {
      }
    }
    const searchSuggestionsState = atom({
      key: `searchSuggestions_${contextID}`,
      default: suggestions,
      effects: [
        ({onSet}) => {
          if (suggestionID) {
            onSet((newValue, _, isReset) => {
              isReset
                ? window.localStorage.removeItem(suggestionID)
                : window.localStorage.setItem(suggestionID, JSON.stringify(newValue))
            })
          }
        }
      ]
    })

    const requiredState = atom({
      key: `required_${contextID}`,
      default: {
        exclude: resource === 'entries' ? ['quantities', 'sections', 'files'] : undefined
      }
    })

    // Stores the available filters. Initialized with static filters, can be
    // modified later.
    const filtersDataState = atom({
      key: `filtersData_${contextID}`,
      default: initialFilterData
    })

    const filterNamesState = selector({
      key: `filters${contextID}`,
      get: ({get}) => new Set(Object.keys(get(filtersDataState)))
    })

    const queryFamily = atomFamily({
      key: `queryFamily_${contextID}`,
      default: (name) => initialQuery[name]
    })

    const dynamicQueryModesFamily = atomFamily({
      key: `dynamicQueryModesFamily_${contextID}`,
      default: (name) => undefined
    })

    // Used to get/set the state of all filters at once
    const filtersState = selector({
      key: `filtersState_${contextID}`,
      get: ({get}) => {
        const query = {}
        for (const key of get(filterNamesState)) {
          const filter = get(queryFamily(key))
          query[key] = filter
        }
        return query
      },
      set: ({get, set}, values) => {
        for (const key of get(filterNamesState)) {
          set(queryFamily(key), values?.[key])
        }
      }
    })

    // Stores a mapping between filter abbreviations and their full names.
    const filterMapsState = selector({
      key: `filterAbbreviations_${contextID}`,
      get: ({get}) => {
        const filterData = get(filtersDataState)
        return getAbbreviations(filterData)
      }
    })

    const paginationState = atom({
      key: `pagination_${contextID}`,
      default: initialPagination
    })

    const guiLocked = parseQueries(initialFiltersLocked, initialFilterData, initialFilterAbbreviations.filterFullnames)
    const lockedFamily = atomFamily({
      key: `lockedFamily_${contextID}`,
      default: (name) => guiLocked?.[name]
    })

    // Used to set the locked state of several filters at once
    const lockedState = selector({
      key: `lockedState_${contextID}`,
      get: ({get}) => {
        const locks = {}
        for (const key of get(filterNamesState)) {
          const filter = get(lockedFamily(key))
          if (!isNil(filter)) locks[key] = filter
        }
        return locks
      }
    })

    /**
     * A Recoil.js selector that aggregates all the currently set filters into a
     * single query object used by the API.
     */
    const queryState = selector({
      key: `query_${contextID}`,
      get: ({get}) => {
        const query = {}
        for (const key of get(filterNamesState)) {
          const filter = get(queryFamily(key))
          if (filter !== undefined) {
            query[key] = filter
          }
        }
        return query
      },
      set: ({ set, get }, data) => {
        for (const filter of get(filterNamesState)) {
          set(queryFamily(filter), undefined)
        }
        if (data) {
          for (const [key, value] of Object.entries(data)) {
            set(queryFamily(key), value)
          }
        }
      }
    })

    /**
     * A Recoil.js selector that return the dynamic query modes state
     */
    const dynamicQueryModesState = selector({
      key: `dynamicQueryModes_${contextID}`,
      get: ({get}) => {
        const query = {}
        for (const key of get(filterNamesState)) {
          const filter = get(dynamicQueryModesFamily(key))
          if (filter !== undefined) {
            query[key] = filter
          }
        }
        return query
      },
      set: ({ set, get }, data) => {
        for (const filter of get(filterNamesState)) {
          set(dynamicQueryModesFamily(filter), undefined)
        }
        if (data) {
          for (const [key, value] of Object.entries(data)) {
            set(dynamicQueryModesFamily(key), value)
          }
        }
      }
    })

    function useDynamicQueryModes() {
      return useRecoilValue(dynamicQueryModesState)
    }

    const widgetFamily = atomFamily({
      key: `widgetFamily_${contextID}`,
      default: (name) => dashboard?.widgets?.[name]
    })

    const widgetIds = atom({
      key: `widgetIds_${contextID}`,
      default: [...Object.keys(dashboard?.widgets || {})]
    })

    // Used to get/set the widgets configuration
    const widgetsState = selector({
      key: `widgetsState_${contextID}`,
      set: ({set, reset}, value) => {
        if (value) {
          const ids = Object.keys(value)
          Object.entries(value).forEach(([key, val]) => set(widgetFamily(key), val))
          set(widgetIds, ids)
        } else {
          reset(widgetIds)
          reset(widgetFamily)
        }
      },
      get: ({get}) => {
        const stats = {}
        const ids = get(widgetIds)
        for (const id of ids) {
          const stat = get(widgetFamily(id))
          if (stat) stats[id] = stat
        }
        return stats
      }
    })

    function useQuery(name) {
      return useRecoilValue(queryState)
    }

    /**
     * This hook will expose a function for reading if widgets are shown
     * on the search page. Use this hook if you intend to only view the value and
     * are not interested in setting the value.
     *
     * @param {string} name Name of the filter.
     * @returns Whether the widgets are shown.
     */
    function useWidgetValue(name) {
      return useRecoilValue(widgetFamily(name))
    }

    /**
     * This hook will expose a function for setting if widget are shown
     * on the search page. Use this hook if you intend to only set the value and are
     * not interested in reading it.
     *
     * @param {string} name Name of the quantity to set.
     * @returns Function for setting the value
     */
    function useSetWidget(name) {
      return useSetRecoilState(widgetFamily(name))
    }

    /**
     * This hook will expose a function for getting and setting whether the
     * widgets are shown. Use this hook if you intend to both read and write the
     * filter value.
     *
     * @param {string} name Name of the filter.
     * @returns Array containing the value and a function for setting it.
     */
    function useWidgetState(name) {
      return useRecoilState(widgetFamily(name))
    }

    /**
     * This hook will expose a function for reading a list of anchored quantities.
     *
     * @returns A list containing the anchored quantity names.
     */
    function useWidgetsValue() {
      return useRecoilValue(widgetsState)
    }

    /**
     * This hook will expose a function for reading and writing the object
     * containing all the current widgets.
     *
     * @returns Array containing the value and a function for setting it.
     */
    function useWidgetsState() {
      return useRecoilState(widgetsState)
    }

    /**
     * Returns a function that can be used to add a new widget.
     *
     * @returns Function for setting a new widget value.
     */
    const useAddWidget = () => {
      const add = useRecoilCallback(({set}) => (id, value) => {
        set(widgetFamily(id), value)
        set(widgetIds, (old) => [...old, id])
      }, [])
      return add
    }

    /**
     * Returns a function that can be used to remove a widget.
     *
     * @returns Function for removing a widget by id.
     */
    const useRemoveWidget = () => {
      const remove = useRecoilCallback(({set}) => (id) => {
        set(widgetFamily(id), undefined)
        set(widgetIds, (old) => [...old].filter(item => item !== id))
      }, [])
      return remove
    }

    /**
     * This hook will expose a function for resetting widgets to the initial
     * configuration.
     *
     * @param {string} name Name of the quantity to set.
     * @returns Function for setting the value
     */
    const useResetWidgets = () => {
      const ids = useRecoilValue(widgetIds)
      const reset = useRecoilCallback(({reset}) => () => {
        for (const id of ids) {
          reset(widgetFamily(id))
        }
        reset(widgetIds)
      }, [])
      return reset
    }

    const aggsFamilyRaw = atomFamily({
      key: `aggsFamilyRaw_${contextID}`,
      default: (name) => initialAggs[name]
    })

    const aggKeys = atom({
      key: `aggKeys_${contextID}`,
      default: []
    })

    const aggsFamily = selectorFamily({
      key: `aggsFamily_${contextID}`,
      get: (id) => ({ get }) => {
        return get(aggsFamilyRaw(id))
      },
      set: (id) => ({set}, value) => {
        set(aggsFamilyRaw(id), value)
        set(aggKeys, prev => [...prev, id])
      }
    })

    /**
     * A Recoil.js selector that aggregates all the currently set aggregation
     * requests into a single query object used by the API.
     */
    const aggsState = selector({
      key: `aggs_${contextID}`,
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
      key: `aggsResponseFamily_${contextID}`,
      default: undefined
    })

    // Recoil.js selector for updating the aggs response in one go.
    const aggsResponseState = selector({
      key: `aggsResponse_${contextID}`,
      get: ({get}) => {
        return undefined
      },
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
      const dynamicQueryModeSetter = useSetRecoilState(dynamicQueryModesFamily(section || name))

      return useCallback((value, config = undefined) => {
        updatedFilters.current.add(name)
        dynamicQueryModeSetter(config?.queryMode)
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
      }, [name, dynamicQueryModeSetter, section, setter, subname])
    }

    /**
     * This hook will expose a function for getting and setting filter values. Use
     * this hook if you intend to both read and write the filter value.
     *
     * @param {string} name Name of the filter.
     * @returns Array containing the filter value and setter function for it.
     */
    const useFilterState = (name) => {
      const value = useFilterValue(name)
      const setter = useSetFilter(name)
      return useMemo(() => [value, setter], [value, setter])
    }

    /**
     * Hook that returns a function that returns the current query as a query string
     *
     * @returns {function} A function that updates the query string to reflect the
     * current search page state.
     */
    const useQueryString = () => {
      const query = useRecoilValue(queryState)
      const filtersData = useRecoilValue(filtersDataState)
      const {filterAbbreviations} = useRecoilValue(filterMapsState)

      return useCallback(() => {
        return convertQueryGUIToQS(query, filtersData, filterAbbreviations)
      }, [query, filtersData, filterAbbreviations])
    }

    /**
     * Returns a function that can be called to reset all current filters.
     *
     * @returns Function for resetting all filters.
     */
    const useResetFilters = () => {
      const filterNames = useRecoilValue(filterNamesState)
      const reset = useRecoilCallback(({set}) => () => {
        for (const filter of filterNames) {
          set(queryFamily(filter), undefined)
          set(dynamicQueryModesFamily(filter), undefined)
        }
      }, [])
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
        const id = `locked_selector_${contextID}_${indexLocked.current}`
        indexLocked.current += 1
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
    const useFilters = (names) => {
      // We dynamically create a Recoil.js selector that is subscribed to the
      // filters specified in the input. This way only the specified filters will
      // cause a render. Recoil.js requires that each selector/atom has an unique
      // id. Because this hook can be called dynamically, we simply generate the ID
      // sequentially.
      const filterState = useMemo(() => {
        const id = `dynamic_selector_${contextID}_${indexFilters.current}`
        indexFilters.current += 1
        return selector({
          key: id,
          get: ({get}) => {
            const query = {}
            for (const key of names) {
              const filter = get(queryFamily(key))
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
      const filtersData = useRecoilValue(filtersDataState)

      // Whenever the aggregation requirements change, create the final
      // aggregation config and set it in the search context: this will then
      // trigger any required API calls that also return the aggregation
      // response that is returned by this hook.
      useEffect(() => {
        const defaults = filtersData[name]?.aggs?.[config?.type]
        const finalConfig = {
          update: update,
          ...defaults || {},
          ...config
        }
        setAgg(finalConfig)
      }, [name, update, setAgg, config, filtersData])

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

    /**
     * Used to update a single filter value with a dynamic name. Supports
     * callback and optionally also setting the queryMode.
     */
    const useUpdateFilter = () => {
      return useRecoilCallback(({set}) => ([key, value, queryMode]) => {
        set(queryFamily(key), value)
        if (queryMode) set(dynamicQueryModesFamily(key), queryMode)
      }, [])
    }

    /**
     * Used to retrieve the currently set search suggestions.
     *
     * @returns Array containing the value and a function for setting it.
     */
    function useSearchSuggestionsValue() {
      return useRecoilValue(searchSuggestionsState)
    }

    /**
     * Used to set the currently set search suggestions.
     */
    function useSetSearchSuggestions() {
      return useSetRecoilState(searchSuggestionsState)
    }

    return {
      useFilterLocked,
      useFiltersLocked,
      useFiltersLockedState,
      useFilterValue,
      useSetFilter,
      useFilterState,
      useFilters,
      useResetFilters,
      useQueryString,
      useDynamicQueryModes,
      aggsFamily,
      useWidgetValue,
      useSetWidget,
      useWidgetState,
      useWidgetsValue,
      useWidgetsState,
      useAddWidget,
      useRemoveWidget,
      useResetWidgets,
      useAgg,
      useSetFilters,
      useUpdateFilter,
      useQuery,
      useSearchSuggestionsValue,
      useSetSearchSuggestions,
      useFilterNames: () => useRecoilValue(filterNamesState),
      useFiltersData: () => useRecoilValue(filtersDataState),
      useFilterMaps: () => useRecoilValue(filterMapsState),
      useSetFiltersData: () => useSetRecoilState(filtersDataState),
      useAggs: () => useRecoilValue(aggsState),
      useSetAggsResponse: () => useSetRecoilState(aggsResponseState),
      usePagination: () => useRecoilValue(paginationState),
      useSetPagination: () => useSetRecoilState(paginationState),
      usePaginationState: () => useRecoilState(paginationState),
      useIsMenuOpen: () => useRecoilValue(isMenuOpenState),
      useSetIsMenuOpen: () => useSetRecoilState(isMenuOpenState),
      useIsCollapsed: () => useRecoilValue(isCollapsedState),
      useSetIsCollapsed: () => useSetRecoilState(isCollapsedState),
      useIsStatisticsEnabled: () => useRecoilValue(isStatisticsEnabledState),
      useSetIsStatisticsEnabled: () => useSetRecoilState(isStatisticsEnabledState),
      useApiQuery: () => useRecoilValue(apiQueryState),
      useSetApiQuery: () => useSetRecoilState(apiQueryState),
      useApiData: () => useRecoilValue(apiDataState),
      useSetApiData: () => useSetRecoilState(apiDataState),
      useResults,
      useSetResults: () => useSetRecoilState(resultsState),
      useResultsUsed: () => useRecoilValue(resultsUsedState),
      useSetResultsUsed: () => useSetRecoilState(resultsUsedState),
      useRequired: () => useRecoilValue(requiredState)
  }
  }, [
    resource,
    contextID,
    initialQuery,
    initialFilterData,
    initialFilterAbbreviations,
    initialFiltersLocked,
    initialPagination,
    dashboard,
    initialAggs,
    id,
    suggestionHistorySize
  ])

  const setResults = useSetResults()
  const setApiData = useSetApiData()
  const setApiQuery = useSetApiQuery()
  const filters = useFilterNames()
  const filtersData = useFiltersData()
  const { filterAbbreviations, filterFullnames } = useFilterMaps()
  const [pagination, setPagination] = usePaginationState()
  const updateAggsResponse = useSetAggsResponse()
  const aggs = useAggs()
  const query = useQuery()
  const dynamicQueryModes = useDynamicQueryModes()
  const filtersLocked = useFiltersLocked()
  const required = useRequired()
  const resultsUsed = useResultsUsed()

  /**
   * This function is used to sync up API calls so that they update the search
   * context state in the same order as they were originally issued.
   *
   * As we cannot guarantee the order in which the API calls finish, we push all
   * calls into a queue. The queue makes sure that API calls get resolved in the
   * original order not matter how long the actual call takes.
   */
  const resolve = useCallback(prop => {
    const {response, timestamp, queryChanged, paginationChanged, search, aggsToUpdate, resource, callbackAgg, callbackHits} = prop
    const data = response.response
    const next = apiQueue.current[0]
    if (next !== timestamp) {
      apiMap.current[timestamp] = prop
      return
    }
    // Update the aggregations if new aggregation data is received. The old
    // aggregation data is preserved and new information is updated.
    if (!isEmpty(data.aggregations)) {
      const newAggs = convertAggAPIToGUI(data.aggregations, aggsToUpdate, resource, filtersData)
      callbackAgg && callbackAgg(newAggs, undefined, true)
    } else {
      callbackAgg && callbackAgg(undefined, undefined, false)
    }
    // Update the query results if new data is received.
    if (queryChanged || paginationChanged) {
      paginationResponse.current = data.pagination
      const newResults = {
        response: response,
        pagination: combinePagination(search.pagination, data.pagination),
        setPagination: setPagination
      }
      callbackHits && callbackHits(newResults, undefined, true, search)
    } else {
      callbackHits && callbackHits(undefined, undefined, false, search)
    }
    // Remove this query from queue and see if next can be resolved.
    apiQueue.current.shift()
    const nextTimestamp = apiQueue.current[0]
    const nextResolve = apiMap.current[nextTimestamp]
    if (nextResolve) {
      resolve(nextResolve)
    }
  }, [setPagination, filtersData])

  /**
   * Function that preprocesses API call requests and finally performs the
   * actual API call.
   *
   * All of the heavier pre-processing, checking, etc. should be done in this
   * function, as it is the final one that gets called after the debounce
   * interval.
   *
   * @param {bool} global Whether the api call is made in the global context
   * that introduces a series of optimizations by keeping track of the previous
   * queries. Should be false when making API calls that should not affect the
   * global context, e.g. when just getting hits with a different pagination etc.
   */
  const apiCall = useCallback((query, aggs, pagination, required, resultsUsed, queryChanged, paginationChanged, updateAggs, callbackAgg = undefined, callbackHits = undefined, global = true) => {
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

    // The locked filters are applied as a parallel AND query. This is the only
    // way to consistently apply them. If we mix them inside 'regular' filters,
    // they can be accidentally overwritten with an OR statement.
    const customQuery = convertQueryGUIToAPI(apiQuery, resource, filtersData, dynamicQueryModes)
    const lockedQuery = convertQueryGUIToAPI(filtersLocked, resource, filtersData, dynamicQueryModes)

    let finalQuery = customQuery
    if (!isEmpty(lockedQuery)) {
      finalQuery = {
        and: [finalQuery, lockedQuery]
      }
    }

    // TODO: This is just a temp hack to somehow allow filter queries with ands/ors
    if (query.custom_quantities) {
      finalQuery = {
        and: [finalQuery, query.custom_quantities]
      }
    }

    const search = {
      owner: filtersLocked?.visibility || apiQuery.visibility,
      query: finalQuery,
      aggregations: convertAggGUIToAPI(
        aggs,
        resource,
        filtersData
      ),
      pagination: convertPaginationGUItoAPI(pagination, filtersData),
      required: convertRequiredGUItoAPI(required, filtersData)
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

    if (global) {
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
    }

    const timestamp = Date.now()
    apiQueue.current.push(timestamp)
    setApiQuery(search?.query)
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
          callbackAgg,
          callbackHits
        })
      })
      .catch((error) => {
        raiseError(error)
        callbackAgg && callbackAgg(undefined, error, true)
        callbackHits && callbackHits(undefined, error, true, undefined)
      })
  }, [filtersData, filterDefaults, filtersLocked, resource, api, raiseError, resolve, dynamicQueryModes, setApiQuery])

  // This is a debounced version of apiCall.
  const apiCallDebounced = useMemo(() => debounce(apiCall, debounceTime), [apiCall])

  /**
   * Intermediate function that should primarily be used when trying to perform
   * an API call. Ensures that ensures that:
   * - Calls are debounced when necessary
   * - API calls are made only if necessary
   */
  const apiCallInterMediate = useCallback((query, aggs, pagination, required, resultsUsed = true, callbackAgg = undefined, callbackHits = undefined, forceUpdate = false, global = true) => {
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
      updatedFilters.current,
      filtersData
    )

    // If results are needed and query and pagination have not changed and
    // aggregations do not need to be updated, no update is necessary. The API
    // calls is made immediately when requesting the first set of results, when
    // the pagination changes or when only aggregations need to be updated.
    // Otherwise it is debounced.
    if ((resultsUsed && (paginationChanged || queryChanged)) || updateAggs) {
      if (firstLoad.current || paginationChanged || !queryChanged) {
        apiCall(query, reducedAggs, pagination, required, resultsUsed, queryChanged, paginationChanged, updateAggs, callbackAgg, callbackHits, global)
      } else {
        apiCallDebounced(query, reducedAggs, pagination, required, resultsUsed, queryChanged, paginationChanged, updateAggs, callbackAgg, callbackHits, global)
      }
    } else {
      callbackAgg && callbackAgg(undefined, undefined, false)
      callbackHits && callbackHits(undefined, undefined, false, undefined)
    }
  }, [apiCall, apiCallDebounced, filtersData])

  // When query, aggregation or pagination changes, update the search context
  useEffect(() => {
    apiCallInterMediate(
      query,
      aggs,
      pagination,
      required,
      resultsUsed,
      (data, error, update) => (update && !error) && updateAggsResponse(data),
      (results, error, update, search) => {
        if (error || !update) return
        setApiData(results.response)
        const data = results.response.response
        const extend = search?.pagination?.page_after_value
        setResults(old => {
          const newResults = old ? {...old} : {}
          extend ? newResults.data = [...newResults.data, ...data.data] : newResults.data = data.data
          newResults.pagination = combinePagination(search.pagination, data.pagination)
          newResults.setPagination = setPagination
          return newResults
        })
      },
      false,
      true
    )
  }, [query, aggs, pagination, required, resultsUsed, apiCallInterMediate, updateAggsResponse, setPagination, setResults, setApiData])

  // The context contains a set of functions that can be used to hook into
  // different data.
  const values = useMemo(() => {
    // Hook for refreshing the results.
    const useRefresh = () => {
      const query = useQuery()
      const aggs = useAggs()
      const pagination = usePagination()
      const required = useRequired()
      const resultsUsed = useResultsUsed()

      const refresh = useCallback(() => {
        apiCallInterMediate(
          query,
          aggs,
          pagination,
          required,
          resultsUsed,
          (data, error, update) => (update && !error) && updateAggsResponse(data),
          (results, error, update, search) => {
            if (error || !update) return
            setApiData(results.response)
            const data = results.response.response
            const extend = search?.pagination?.page_after_value
            setResults(old => {
              const newResults = old ? {...old} : {}
              extend ? newResults.data = [...newResults.data, ...data.data] : newResults.data = data.data
              newResults.pagination = combinePagination(search.pagination, data.pagination)
              newResults.setPagination = setPagination
              return newResults
            })
          },
          true,
          true
        )
      }, [aggs, pagination, required, query, resultsUsed])
      return refresh
    }

    // Hook for imperatively requesting aggregation data. By using this hook you
    // can track the state of individual calls and perform callbacks.
    const useAggCall = (name, id) => {
      const key = useMemo(() => `${name}:${id}`, [name, id])
      const query = useQuery()
      const pagination = usePagination()
      const required = useRequired()
      const resultsUsed = useResultsUsed()
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
          query,
          aggs,
          pagination,
          required,
          resultsUsed,
          (response, error, update) => {
            if (update && !error) {
              updateAggsResponse(response)
            }
            callback(response && response[key], error, update)
          },
          undefined,
          false,
          true
        )

        // We also need to update aggregation request state, otherwise the
        // subsequent calls will not be able to know what was done by this call.
        // To do this without triggering another API call, we disable API updates
        // for one cycle.
        disableUpdate.current = true
        aggs[key].update = false
        setAgg(aggs[key])
      }, [key, query, pagination, required, resultsUsed, setAgg])

      return aggCall
    }

    /**
     * Hook for fetching an independent set of hits based on the currently set
     * filters.
     *
     * @param {string} id Unique identified for this set of hits
     * @param {object} required The required field of the API call
     * @param {object} pagination The API pagination
     *
     * @returns {Array} Array of hit objects.
     */
    const useHits = (id = 'default', required, pagination, callback) => {
      const [results, setResults] = useState()
      const query = useQuery()
      const apiCallDebounced = useMemo(() => debounce(apiCall, debounceTime), [])

      // Update hits when query, includes, excludes or pagination change.
      useEffect(() => {
        apiCallDebounced(
          query,
          {},
          convertPaginationGUItoAPI(pagination, filtersData),
          convertRequiredGUItoAPI(required, filtersData),
          true,
          true,
          true,
          false,
          undefined,
          (results, error, update) => {
            if (!error && update) {
              const data = results.response.response
              update && !error && setResults(data?.data)
            }
            callback(results, error, update)
          },
          false
        )
      }, [query, required, pagination, apiCallDebounced, callback])

      return results
    }

    /**
     * Hook that can be used to fetch search suggestions stored in the local
     * device.
     * */
    const useSearchSuggestions = (input, size = 5) => {
      const suggestions = useSearchSuggestionsValue()
      const suggestionFilter = useMemo(() => {
        const sugg = getSuggestions(suggestions, 0, suggestion => suggestion.input)
        return sugg.filter
      }, [suggestions])

      return input.trim() === ''
        ? suggestions.slice(0, size)
        : suggestionFilter(input).slice(0, size)
    }

    /**
     * Hook for pushing new suggestion values to the search context. New values
     * are pushed to the top of the list, and if the suggestion already exists,
     * it is simply pulled to the top.
     * */
    const usePushSearchSuggestion = () => {
      const setSuggestions = useSetSearchSuggestions()
      return useCallback((suggestion) => {
        setSuggestions(old => {
          let newSuggestions = old.filter(x => x.key !== suggestion.key)
          const newSuggestion = new SearchSuggestion(suggestion)
          newSuggestion.history = true
          newSuggestions = [newSuggestion, ...old].slice(0, suggestionHistorySize)
          return newSuggestions
        })
      }, [setSuggestions])
    }

    /**
     * Hook for deleting specific suggestion values from the search context
     * using the suggestion key.
     * */
    const useRemoveSearchSuggestion = () => {
      const setSuggestions = useSetSearchSuggestions()
      return useCallback((key) => {
        setSuggestions(old => {
          const newSuggestions = [...old].filter(x => x.key !== key)
          return newSuggestions
        })
      }, [setSuggestions])
    }

    /**
     * Hook that returns a function for parsing a single filter key, value pair
     * into a form that is supported by the GUI. The parsing depends on the
     * SearchContext as filters may define a custom deserialization function and
     * only certain filters may be available.
     * */
    const useParseQuery = () => {
      return useCallback(
        (key, value, path, multiple) => parseQuery(key, value, filtersData, units, path, multiple),
        []
      )
    }
    /**
     * Hook that returns a function for parsing an object containing filter
     * values into a form that is supported by the GUI. The parsing depends on
     * the SearchContext as filters may define a custom deserialization function
     * and only certain filters may be available.
     * */
    const useParseQueries = () => {
      return useCallback(
        (values) => parseQueries(values, filtersData, filterFullnames),
        []
      )
    }

    return {
      resource,
      columns,
      rows,
      filterMenus,
      filters,
      filterData: filtersData,
      filterGroups: initialFilterGroups,
      filterFullnames,
      filterAbbreviations,
      searchSyntaxes: initialSearchSyntaxes,
      useIsMenuOpen,
      useSetIsMenuOpen,
      useIsCollapsed,
      useSetIsCollapsed,
      useIsStatisticsEnabled,
      useSetIsStatisticsEnabled,
      useApiData,
      useApiQuery,
      useQuery,
      useSetPagination,
      useParseQuery,
      useParseQueries,
      useSearchSuggestions,
      usePushSearchSuggestion,
      useRemoveSearchSuggestion,
      useFilterValue,
      useSetFilter,
      useFilterState,
      useFilters,
      useResetFilters,
      useRefresh,
      useFilterLocked,
      useFiltersLocked,
      useFiltersLockedState,
      useWidgetValue,
      useSetWidget,
      useWidgetState,
      useWidgetsValue,
      useWidgetsState,
      useAddWidget,
      useRemoveWidget,
      useResetWidgets,
      useQueryString,
      useResults,
      useHits,
      useAgg,
      useAggCall,
      useSetFilters,
      useUpdateFilter
    }
  }, [
    resource,
    rows,
    columns,
    filterMenus,
    filters,
    filtersData,
    initialFilterGroups,
    filterFullnames,
    filterAbbreviations,
    initialSearchSyntaxes,
    useIsMenuOpen,
    useSetIsMenuOpen,
    useIsCollapsed,
    useSetIsCollapsed,
    useIsStatisticsEnabled,
    useSetIsStatisticsEnabled,
    useFilterValue,
    useSetFilter,
    useFilterState,
    useFilters,
    useResetFilters,
    useFilterLocked,
    useFiltersLocked,
    useFiltersLockedState,
    useWidgetValue,
    useSetWidget,
    useWidgetState,
    useWidgetsValue,
    useWidgetsState,
    useAddWidget,
    useRemoveWidget,
    useResetWidgets,
    useQueryString,
    usePagination,
    useSetPagination,
    useResults,
    useResultsUsed,
    useRequired,
    useApiQuery,
    useApiData,
    useAgg,
    useAggs,
    useQuery,
    useSetFilters,
    useUpdateFilter,
    apiCallInterMediate,
    apiCall,
    aggsFamily,
    updateAggsResponse,
    setPagination,
    setResults,
    setApiData,
    units,
    suggestionHistorySize,
    useSearchSuggestionsValue,
    useSetSearchSuggestions
  ])

  return <searchContext.Provider value={values}>
    {children}
  </searchContext.Provider>
})

SearchContextRaw.propTypes = {
  resource: PropTypes.string,
  initialFiltersLocked: PropTypes.object,
  initialColumns: PropTypes.object,
  initialRows: PropTypes.object,
  initialFilterMenus: PropTypes.object,
  initialPagination: PropTypes.object,
  initialDashboard: PropTypes.object,
  initialFilterData: PropTypes.object, // Determines which filters are available
  initialFilterGroups: PropTypes.object, // Maps filter groups to a set of filter names
  initialFilterValues: PropTypes.object, // Here one can provide default filter values
  initialSearchSyntaxes: PropTypes.object, // Determines which syntaxes are supported
  children: PropTypes.node,
  id: PropTypes.string,
  suggestionHistorySize: PropTypes.number
}

SearchContextRaw.defaultProps = {
  suggestionHistorySize: 20
}

export const SearchContext = compose(withQueryString, withFilters)(SearchContextRaw)

/**
 * Hook for accessing the current SearchContext.
 */
export function useSearchContext() {
  return useContext(searchContext)
}

/**
 * Parses a single filter value into a form that is supported by the GUI. This includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 * - Custom serialization defined for the filter will be used if defined.
 *
 * @param {string} key Name of the filter or an operator name.
 * @param {any} value Value of the filter or operator
 * @param {object} filtersData All of the filters that are available
 * @param {object} units Unit system for unit conversion
 * @param {string} path The full path for the filter without an operator name.
 *
 * @returns {any} The filter value in a format that is suitable for the GUI.
 */
function parseQuery(key, value, filtersData, units = undefined, path = undefined, multiple = undefined) {
  let newValue
  const fullPath = path ? `${path}.${key}` : key
  if (isPlainObject(value)) {
    newValue = {}
    for (const [keyInner, valueInner] of Object.entries(value)) {
      const valueConverted = parseQuery(keyInner, valueInner, filtersData, units, fullPath)
      if (!isNil(valueConverted)) {
        newValue[keyInner] = valueConverted
      }
    }
  } else {
    // If the key is an operator, the filter name is read from the path.
    const opKeys = new Set(['lte', 'lt', 'gte', 'gt'])
    const isOperator = opKeys.has(key)
    const filterPath = isOperator ? path : fullPath
    multiple = multiple ?? filtersData[filterPath].multiple
    const deserializer = filtersData[filterPath].deserializer
    if (isArray(value) || isSet(value)) {
      newValue = new Set(value.map((v) => deserializer(v, units)))
    } else {
      newValue = deserializer(value, units)
      if (!isOperator && multiple) {
        newValue = new Set([newValue])
      }
    }
  }
  return newValue
}

/**
 * Cleans an entire query object into a form that is supported by the GUI. This
 * includes:
 * - Arrays are are transformed into Sets
 * - If multiple values are supported, scalar values are stored inside sets.
 * - Numerical values with units are transformed into Quantities.
 *
 * @param {object} query Query object to transform.
 * @param {object} filtersData All of the filters that are available
 * @param {object} units The desired unit system used when reading quantities.
 *
 * @returns {any} The filter object in a format that is suitable for the GUI.
 */
function parseQueries(query, filtersData, filterFullnames, units = undefined) {
  const newQuery = {}
  if (query) {
    for (const [key, value] of Object.entries(query)) {
      const newKey = filterFullnames[key] || key
      const valueGUI = parseQuery(newKey, value, filtersData, units)
      newQuery[newKey] = valueGUI
    }
  }
  return newQuery
}

/**
 * Converts a query object into a valid query string.
 *
 * @param {object} query Query object representing the currently active
 * filters in a search context.
 * @param {object} filtersData All of the filters that are available
 * @returns URL querystring, not encoded to improve readability.
 */
function convertQueryGUIToQS(query, filtersData, filterAbbreviations) {
  // Used to recursively convert the query into a serializable format.
  function convert(key, value, path) {
    // If the key is an operator, the filter name is read from the path.
    const opKeys = new Set(['lte', 'lt', 'gte', 'gt'])
    const fullPath = path ? `${path}.${key}` : key
    const filterPath = opKeys.has(key) ? path : fullPath
    const filterData = filtersData[filterPath]
    let newValue
    if (isPlainObject(value) && !filterData.customSerialization) {
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
        queryStringQuery[key] = valueConverted
      }
    }
  }

  return qs.stringify(queryStringQuery, {indices: false, encode: false})
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
 * @param {object} filtersData All of the filters that are available
 * @param {object} queryModes Query modes for filters
 *
 * @returns {object} A copy of the object with certain items cleaned into a
 * format that is supported by the API.
 */
function convertQueryGUIToAPI(query, resource, filtersData, queryModes) {
  const queryCustomized = {}
  if (!query) {
    return undefined
  }

  // Perform custom transformations. Note that these calls may add new keys to
  // the query.
  function customize(key, value, parent, subKey = undefined) {
    const data = filtersData[key]

    // Global filters are not serialized into the API call.
    if (data?.global) return

    // Sections need to be recursively handled. Notice that we can't directly
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
    if (!isNil(v)) {
      customize(k, v, queryCustomized)
    }
  }

  // Create the API-compatible keys and values.
  const queryNormalized = {}
  for (const [k, v] of Object.entries(queryCustomized)) {
    const newValue = convertQuerySingleGUIToAPI(v)
    const {quantity: filterName, op: queryMode} = parseOperator(k)
    let finalKey = filtersData[filterName]?.requestQuantity || filterName
    finalKey = resource === 'materials' ? getFilterMaterialPath(finalKey) : finalKey
    let finalQueryMode = queryMode || queryModes?.[k]
    if (isNil(finalQueryMode) && isArray(newValue)) {
      finalQueryMode = filtersData[k]?.queryMode
    }
    finalKey = finalQueryMode ? `${finalKey}:${finalQueryMode}` : finalKey
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
          const {quantity: newKey, op: queryMode} = parseOperator(k)
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
 * - Empty containers are set to undefined
 *
 * @param {string} key Filter name
 * @param {any} value Filter value
 * @param {object} filtersData All of the filters that are available
 * @param {string} queryMode Determines the queryMode
 *
 * @returns {any} The filter value in a format that is suitable for the API.
 */
function convertQuerySingleGUIToAPI(value) {
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
      const apiValue = convertQuerySingleGUIToAPI(valueInner)
      if (!isNil(apiValue)) {
        newValue[keyInner] = apiValue
      }
    }
    if (isEmpty(newValue)) {
      newValue = undefined
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
 * @param {object} aggs The aggregation data as constructed by the GUI.
 * @param {object} updatedFilters Set of filters that were updated together with
 * this call.
 * @param {string} resource The resource we are looking at: entries or materials.
 * @param {object} filtersData All of the filters that are available
 *
 * @returns {object} Aggregation query that is usable by the API.
 */
function convertAggGUIToAPI(aggs, resource, filtersData) {
  const apiAggs = {}
  for (const [key, agg] of Object.entries(aggs)) {
    const filterName = rsplit(key, ':', 1)[0]
    if (agg.update) {
      const exclusive = filtersData[filterName].exclusive
      const type = agg.type
      const apiAgg = apiAggs[key] || {}
      const aggSet = agg.set
      const finalAgg = aggSet ? aggSet(agg) : agg
      let finalQuantity = finalAgg.quantity || filterName

      // If the targeted quantity has a special requestQuantity, it is used in
      // the final API call.
      finalQuantity = filtersData[finalQuantity]?.requestQuantity || finalQuantity

      // The targeted quantity is decided based on the resource as materials
      // search need to target slightly different fields.
      finalQuantity = resource === 'materials'
        ? getFilterMaterialPath(finalQuantity)
        : finalQuantity

      apiAgg[type] = {
        // Exclusive quantities (quantities that have one value per entry) are
        // always fetched with exclude_from_search
        exclude_from_search: exclusive,
        ...finalAgg,
        quantity: finalQuantity
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
 * @param {object} filtersData All of the filters that are available
 *
 * @returns {object} Aggregation result that is usable by the GUI.
 */
function convertAggAPIToGUI(aggs, aggsToUpdate, filtersData) {
  if (isEmpty(aggs)) {
    return {}
  }

  // Perform custom transformations
  const aggsCustomized = {}
  for (const key of aggsToUpdate) {
    const filter_name = rsplit(key, ':', 1)[0]
    const aggConfig = aggs[key]
    if (!isNil(aggConfig)) {
      for (const [type, agg] of Object.entries(aggConfig)) {
        const aggGet = filtersData[filter_name]?.aggs?.[type]?.get
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
 * Used to transform a GUI pagination call result into a form that is usable by the
 * API.
 *
 * @param {object} pagination The aggregation data as returned by the API.
 * @param {object} filtersData All of the filters that are available
 *
 * @returns {object} Pagination object that is usable by the API.
 */
function convertPaginationGUItoAPI(pagination, filtersData) {
  if (!pagination) return {}
  return {
    ...pagination,
    order_by: filtersData[pagination.order_by]?.requestQuantity || pagination.order_by
  }
}

/**
 * Used to transform a GUI pagination call result into a form that is usable by the
 * API.
 *
 * @param {object} pagination The aggregation data as returned by the API.
 * @param {object} filtersData All of the filters that are available
 *
 * @returns {object} Pagination object that is usable by the API.
 */
function convertRequiredGUItoAPI(required, filtersData) {
  if (!required) return {}
  return {
    ...required,
    include: required.include?.map(name => filtersData[name]?.requestQuantity || name),
    exclude: required.exclude?.map(name => filtersData[name]?.requestQuantity || name)
  }
}

/**
 * Goes through the given aggregations and compares then against the old results
 * and the current query to identify and return the aggregations that need to be
 * called.
 *
 * @param {object} aggs The current aggregation configuration.
 * @param {object} oldAggs Reduced aggregation config from latest finished query.
 * @param {bool} queryChanged Whether the query has changed.
 * @param {object} filtersData All of the filters that are available
 *
 * @returns {object} Reduced aggregation config.
 */
function reduceAggs(aggs, oldAggs, queryChanged, updatedFilters, filtersData) {
  // Loop through the different aggregations for each quantity and see if any
  // of them need to be updated.
  const reducedAggs = {}
  let updateAggs = false
  for (const [key, agg] of Object.entries(aggs)) {
    const filter_name = rsplit(key, ':', 1)[0]
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
      ? filtersData[filter_name].exclusive
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

/**
 * Returns the final value that should be shown for the given filter.
 *
 * @param {*} def The quantity definition
 * @param {*} filter The current filter value in the search context
 * @param {*} filterLocked The current locked filter value in the search context
 * @param {*} initialValue Initial value that overrides any default specified in
 *   the FilterRegistry
 *
 * @returns The final value for the filter.
 */
export function getValue(def, filter, filterLocked, initialValue) {
  const defaultValue = isNil(initialValue) ? def.default : initialValue
  const val = (def.global && !isNil(filterLocked)) ? filterLocked : filter
  return isNil(val) ? defaultValue : val
}

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

/**
 * Returns a mapping from full filter names to their abbreviations
 * and vice-versa.
*/
function getAbbreviations(filterData) {
  const abbreviations = {}
  const filterAbbreviations = {}
  const filterFullnames = {}
  const nameAbbreviationPairs = Object.keys(filterData).map(
    fullname => {
      const {path, schema} = parseQuantityName(fullname)
      return [fullname, schema ? path : path.split('.').pop()]
    }
  )
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
  return {filterAbbreviations, filterFullnames}
}

/**
 * Function that is used to automatically convert paths that target the correct
 * property in the materials index.
 * @param {string} path Path in the entry index
 * @returns The correct path in the materials index
 */
function getFilterMaterialPath(path) {
  const prefix = 'results.material.'
  let materialPath
  if (path.startsWith(prefix)) {
    materialPath = path.substring(prefix.length)
  } else {
    materialPath = `entries.${path}`
  }
  return materialPath
}
