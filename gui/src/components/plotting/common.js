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
import { range, isNil } from 'lodash'
import { scalePow } from 'd3-scale'
import {
  format,
  getTime,
  fromUnixTime,
  millisecondsToSeconds,
  eachMinuteOfInterval,
  eachHourOfInterval,
  eachDayOfInterval,
  eachWeekOfInterval,
  eachMonthOfInterval,
  eachYearOfInterval,
  differenceInYears,
  differenceInMonths,
  differenceInWeeks,
  differenceInDays,
  differenceInHours,
  differenceInMinutes,
  differenceInSeconds,
  isWithinInterval,
  eachQuarterOfInterval
} from 'date-fns'
import { scale, add, DType, formatNumber } from '../../utils.js'

// The available scaling options
export const scales = {
  'linear': 1,
  '1/2': 0.5,
  '1/4': 0.25,
  '1/8': 0.125
}

/**
 * Returns a d3-scale object for the given scaling type, domain and range.
 *
 * @param {str} type Type of scaling.
 * @param {array} domain The input domain.
 * @param {array} range The output range.
 * @returns A function that when given a number between [0, 1] will transform
 * it to the given output range with the given type of scaling.
 */
export function getScaler(type, domain = [0, 1], range = [0, 1]) {
  const scale = scales[type]
  if (isNil(scale)) {
    throw Error('Invalid scaling type.')
  }
  const scaler = scalePow()
    .exponent(scale)
    .domain(domain)
    .range(range)

  return scaler
}

/**
 * Returns interval that can be used to split the given range into
 * human-readable, uniformly spaced steps.
 *
 * @param {number} range The range to divide
 * @param {number} steps The targeted number of steps.
 * @param {string} dtype Data type.
 * @param {bool} cap controls whether the final number of intervals can be
 * bigger than the requested amount.
 *
 * @returns The interval as a number
 */
export function getInterval(value, steps, dtype, cap = true) {
  const interval = value / steps
  const degree = Math.pow(10, (Math.round(Math.log10(interval))))
  const multipliers = [0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10]
  const intervals = scale(multipliers, degree).filter(x => {
    const capped = cap ? (x >= interval) : true
    const valid = dtype === DType.Int ? x % 1 === 0 : true
    return capped && valid
  })
  const differences = add(intervals, -interval).map(x => Math.abs(x))
  const indexMin = differences.reduce(argMin, 0)
  return intervals[indexMin]
}

export const argMin = (iMin, x, i, arr) => {
  return x < arr[iMin] ? i : iMin
}

/**
 * Returns human-readable, (mostly) uniformly spaced ticks for the given
 * interval. Takes the data type into account when determinining the ticks their
 * formatting. If not valid ticks can be generated, returns an empty list.
 *
 * @param {number} min The minimum value
 * @param {number} max The maximum value
 * @param {number} n Targeted number of ticks. Notice that this typically does
 * not correspond to the final number of ticks.
 * @param {string} dtype Data type for the values
 * @param {number} decimals Number of decimals to show
 * @param {bool} scientific Whether to use scientific notation
 *
 * @returns Array of tick objects containing value and tick.
 */
export function getTicks(min, max, n, dtype, mode = 'scientific', decimals = 3) {
  if (dtype === DType.Timestamp) {
    const start = fromUnixTime(millisecondsToSeconds(min))
    const end = fromUnixTime(millisecondsToSeconds(max))
    const interval = {start, end}

    // Function for splitting ranges using milliseconds
    const splitMilliseconds = (step) => ({start, end}) => {
      const startRound = Math.ceil(getTime(start) / step)
      const endRound = Math.floor(getTime(end) / step)
      return range(startRound, endRound + 1).map(x => x * step)
    }

    // Config for the available intervals
    const intervals = {
      years: {
        difference: differenceInYears,
        split: eachYearOfInterval,
        format: 'yyyy'
      },
      threemonths: {
        difference: (end, start) => differenceInMonths(end, start) / 3,
        split: eachQuarterOfInterval,
        format: 'MMM yyyy'
      },
      months: {
        difference: differenceInMonths,
        split: eachMonthOfInterval,
        format: 'MMM'
      },
      weeks: {
        difference: differenceInWeeks,
        split: eachWeekOfInterval,
        format: 'MMM d'
      },
      days: {
        difference: differenceInDays,
        split: eachDayOfInterval,
        format: 'MMM d'
      },
      twelvehours: {
        difference: (end, start) => differenceInHours(end, start) / 12,
        split: (interval) => eachHourOfInterval(interval).filter(x => {
          const values = new Set(['12:00', '24:00'])
          return values.has(format(x, 'k:mm'))
        }),
        format: 'k:mm'
      },
      sixhours: {
        difference: (end, start) => differenceInHours(end, start) / 6,
        split: (interval) => eachHourOfInterval(interval).filter(x => {
          const values = new Set(['6:00', '12:00', '18:00', '24:00'])
          return values.has(format(x, 'k:mm'))
        }),
        format: 'k:mm'
      },
      hours: {
        difference: differenceInHours,
        split: eachHourOfInterval,
        format: 'k:mm'
      },
      thirtyminutes: {
        difference: (end, start) => differenceInMinutes(end, start) / 30,
        split: (interval) => eachMinuteOfInterval(interval, {step: 30}),
        format: 'k:mm'
      },
      fifteenminutes: {
        difference: (end, start) => differenceInMinutes(end, start) / 15,
        split: (interval) => eachMinuteOfInterval(interval, {step: 15}),
        format: 'k:mm'
      },
      minutes: {
        difference: differenceInMinutes,
        split: eachMinuteOfInterval,
        format: 'k:mm'
      },
      thirtyseconds: {
        difference: (end, start) => differenceInSeconds(end, start) / 30,
        split: splitMilliseconds(30000),
        format: 'k:mm:ss'
      },
      fifteensecond: {
        difference: (end, start) => differenceInSeconds(end, start) / 15,
        split: splitMilliseconds(15000),
        format: 'k:mm:ss'
      },
      seconds: {
        difference: differenceInSeconds,
        split: splitMilliseconds(1000),
        format: 'k:mm:ss'
      }
    }

    // Calculate the minimum number of ticks that each option may produce. The
    // actual number of ticks cannot be easily calculated without actually
    // producing the ticks, which is not viable when there may be a very large
    // number of them (e.g. milliseconds in a year.)
    const differences = {}
    Object.entries(intervals)
      .forEach(([dur, setup]) => { differences[dur] = setup.difference(end, start) })

    // For each sensible option, calculate the actual ticks and their amount.
    const differencesValid = Object.entries(differences)
      .filter(([dur, diff]) => diff <= n) // Filter out values where the minimum number of ticks is beyond the target.
      .map(([dur, diff]) => [dur, intervals[dur].split(interval).filter(x => isWithinInterval(x, interval))]) // Calculate the actual ticks
      .filter(([dur, ticks]) => ticks.length <= n) // Filter out options with more ticks than requested
      .map(([dur, ticks]) => [dur, Math.abs(ticks.length - n), ticks]) // Calculate abs difference to n

    // If no valid ticks can be determined, return an empty list.
    if (!differencesValid.length) {
      return []
    }

    // Determine the option that has the closest amount of ticks to the targeted
    // amount.
    const closestDuration = differencesValid
      .reduce((prev, curr) => prev[1] < curr[1] ? prev : curr)
    const ticks = closestDuration[2]
    const setup = intervals[closestDuration[0]]

    // Format the results and return both the formatted tick and the final
    // value.
    return ticks.map(x => ({value: x, tick: format(x, setup.format)}))
  } else {
    // Calculate minimum number of ticks for each option
    const tickRange = max - min
    const multipliers = [0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10]
    const degree = Math.pow(10, (Math.round(Math.log10(tickRange))))
    const closestDuration = scale(multipliers, degree)
      .filter(x => dtype === DType.Int ? x % 1 === 0 : true) // Filter out invalid intervals
      .map(x => [x, Math.floor(tickRange / x)]) // Calculate minimum number of ticks
      .filter(([interval, nMin]) => nMin <= n) // Filter out values where the minimum number of ticks is beyond the target.
      .map(([interval, nMin]) => { // Calculate actual ticks
        const startRound = Math.ceil(min / interval)
        const endRound = Math.floor(max / interval)
        const ticks = range(startRound, endRound + 1).map(x => {
          const value = (x * interval)
          return {value, tick: formatNumber(value, dtype, mode, decimals)}
        })
        return ticks
      })
      .filter((ticks) => ticks.length <= n) // Filter out options with more ticks than requested
      .map((ticks) => [ticks, Math.abs(ticks.length - n)]) // Calculate abs difference to n
      .reduce((prev, curr) => prev[1] < curr[1] ? prev : curr) // Select best option
    return closestDuration[0]
  }
}
