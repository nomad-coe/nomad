import React, { useContext, useMemo, useCallback } from 'react'
import PropTypes from 'prop-types'
import { searchContext } from './SearchContext.js'
import searchQuantities from '../../searchQuantities'
import Histogram from '../Histogram.js'

const unprocessedLabel = 'not processed'
const unavailableLabel = 'unavailable'

export default function QuantityHistogram({
  quantity, valueLabels = {}, title, values, numberOfValues, multiple, tooltips = {},
  ...props
}) {
  title = title || quantity
  values = values || (searchQuantities[quantity] && searchQuantities[quantity].statistic_values)
  numberOfValues = numberOfValues || (values && values.length) || (searchQuantities[quantity] && searchQuantities[quantity].statistic_size)
  const {response: {statistics, metric}, query, setQuery} = useContext(searchContext)
  const statisticsData = statistics[quantity]

  const handleItemClicked = useCallback(item => {
    console.log(quantity)
    console.log(query)
    if (multiple) {
      // Add or remove item from query
      let newQuery = query[quantity]
      console.log("Old query: " + newQuery)
      if (newQuery === undefined) {
        newQuery = [item.key]
      } else {
        if (!Array.isArray(newQuery)) {
          newQuery = [newQuery]
        }
        newQuery = new Set(newQuery)
        if (newQuery.has(item.key)) {
          newQuery.delete(item.key)
        } else {
          newQuery.add([item.key])
        }
        newQuery = Array.from(newQuery.values())
      }
      console.log("New query: " + newQuery)
      setQuery({[quantity]: newQuery})
    } else {
      setQuery({[quantity]: (query[quantity] === item.key) ? null : item.key})
    }
  }, [query, setQuery, multiple, quantity])

  const data = useMemo(() => {
    let data
    if (!statistics[quantity]) {
      data = []
    } else if (values) {
      data = values.map(value => ({
        key: value,
        name: valueLabels[value] || value,
        value: statisticsData[value] ? statisticsData[value][metric] : 0,
        tooltip: tooltips[value]
      }))
    } else {
      data = Object.keys(statisticsData)
        .map(value => ({
          key: value,
          name: valueLabels[value] || value,
          value: statisticsData[value][metric]
        }))
      // keep the data sorting, but put unavailable and not processed to the end
      const unavailableIndex = data.findIndex(d => d.name === unavailableLabel)
      const unprocessedIndex = data.findIndex(d => d.name === unprocessedLabel)
      if (unavailableIndex !== -1) {
        data.push(data.splice(unavailableIndex, 1)[0])
      }
      if (unprocessedIndex !== -1) {
        data.push(data.splice(unprocessedIndex, 1)[0])
      }
    }
    return data
  }, [metric, quantity, statistics, statisticsData, valueLabels, values, tooltips])

  return <Histogram
    card data={data}
    numberOfValues={numberOfValues}
    title={title}
    onClick={handleItemClicked}
    selected={query[quantity]}
    multiple={multiple}
    tooltips={!!tooltips}
    {...props}
  />
}
QuantityHistogram.propTypes = {
  /**
   * The name of the search quantity that is displayed in the histogram. This has to
   * match the provided statistics data.
   */
  quantity: PropTypes.string.isRequired,
  /**
   * An optional title for the chart. If no title is given, the quantity is used.
   */
  title: PropTypes.string,
  /**
   * The data. Usually the statistics data send by NOMAD's API.
   */
  data: PropTypes.object,
  /**
   * Optional list of possible values. This is used to sort the data and fill the data
   * with 0-values to keep a persistent appearance, even if no data for that value exists.
   * Otherwise, the values are not sorted.
   */
  values: PropTypes.arrayOf(PropTypes.string),
  /**
   * The maximum number of values. This is used to fix the histograms size. Otherwise,
   * the size is determined by the required space to render the existing values.
   */
  numberOfValues: PropTypes.number,
  /**
   * An optional mapping between values and labels that should be used to render the
   * values.
   */
  valueLabels: PropTypes.object,
  /**
   * An optional mapping between values and their tooltip content.
   */
  tooltips: PropTypes.object,
  /**
   * Whether multiple values can be appended to the same query key.
   */
  multiple: PropTypes.bool
}

QuantityHistogram.defaultProps = {
  multiple: false
}
