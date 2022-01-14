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
import { cloneDeep, merge, isSet, isNil, isString, isNumber } from 'lodash'
import { toUnitSystem, Quantity } from './units'
import { fromUnixTime, format } from 'date-fns'
import { dateFormat } from './config'
import { scale as chromaScale } from 'chroma-js'
import searchQuantities from './searchQuantities.json'

export const isEquivalent = (a, b) => {
  // Create arrays of property names
  var aProps = Object.getOwnPropertyNames(a)
  var bProps = Object.getOwnPropertyNames(b)

  // If number of properties is different,
  // objects are not equivalent
  if (aProps.length !== bProps.length) {
    return false
  }

  for (var i = 0; i < aProps.length; i++) {
    var propName = aProps[i]

    // If values of same property are not equal,
    // objects are not equivalent
    if (a[propName] !== b[propName]) {
      return false
    }
  }

  // If we made it this far, objects
  // are considered equivalent
  return true
}

export const capitalize = (s) => {
  if (typeof s !== 'string') {
    return ''
  }
  return s.charAt(0).toUpperCase() + s.slice(1)
}

/**
 * Used to scale numeric values. Works on n-dimensional arrays and implemented
 * as a relatively simple for loop for performance. If conversion times become
 * an issue, it might be worthwhile to look at vectorization with WebAssembly.
 *
 * @param {*} value The original values.
 * @param {number} addition Number to add.
 *
 * @return {*} A copy of the original data with numbers scaled.
 */
export function scale(value, factor) {
  // Function for recursively scaling the values
  function scaleRecursive(list, newList) {
    const isScalarArray = !Array.isArray(list[0])
    if (isScalarArray) {
      for (let i = 0, size = list.length; i < size; ++i) {
        newList.push(list[i] * factor)
      }
    } else {
      for (let i = 0, size = list.length; i < size; ++i) {
        const iList = []
        newList.push(iList)
        scaleRecursive(list[i], iList)
      }
    }
  }
  // If given a scalar variable, simply try to scale it. If a list is given,
  // perform the scaling recursively.
  const isArray = Array.isArray(value)
  let newValue
  if (!isArray) {
    newValue = value * factor
  } else {
    newValue = []
    scaleRecursive(value, newValue)
  }
  return newValue
}

/**
 * Used to add a single scalar value to an n-dimensional array.
 *
 * @param {*} value The original values.
 * @param {number} addition Number to add.
 *
 * @return {*} A copy of the original data with the addition made for each
 * number.
 */
export function add(value, addition) {
  // Function for recursively adding the value
  function addRecursive(list, newList) {
    const isScalarArray = !Array.isArray(list[0])
    if (isScalarArray) {
      for (let i = 0, size = list.length; i < size; ++i) {
        newList.push(list[i] + addition)
      }
    } else {
      for (let i = 0, size = list.length; i < size; ++i) {
        const iList = []
        newList.push(iList)
        addRecursive(list[i], iList)
      }
    }
  }
  // If given a scalar variable, simply try to make the addition. If a list is
  // given, perform the addition recursively.
  const isArray = Array.isArray(value)
  let newValue
  if (!isArray) {
    newValue = value + addition
  } else {
    newValue = []
    addRecursive(value, newValue)
  }
  return newValue
}

/**
 * Used to calculate the distance between two n-dimensional points,
 *
 * @param {*} a First point
 * @param {*} b Second point
 *
 * @return {*} Euclidean distance between the given two points.
 */
export function distance(a, b) {
  return a
    .map((x, i) => Math.abs(x - b[i]) ** 2) // square the difference
    .reduce((sum, now) => sum + now) ** // sum
    (1 / 2)
}

/**
 * Used to merge two Javascript objects into a new third object by recursively
 * overwriting and extending the target object with properties from the source
 * object.
 *
 * @param {*} target The values to convert
 * @param {*} source Original unit.
 *
 * @return {*} A copy of the original data with units converted.
 */
export function mergeObjects(source, target, copy = false) {
  // First create a deep clone that will be used as the returned object
  let cloned = cloneDeep(target)
  let val = merge(cloned, source)
  return val
}

export function arraysEqual(_arr1, _arr2) {
  if (!Array.isArray(_arr1) || !Array.isArray(_arr2) || _arr1.length !== _arr2.length) {
    return false
  }

  var arr1 = _arr1.concat().sort()
  var arr2 = _arr2.concat().sort()

  for (var i = 0; i < arr1.length; i++) {
    if (arr1[i] !== arr2[i]) { return false }
  }

  return true
}

export function onlyUnique(value, index, self) {
  return self.indexOf(value) === index
}

export function objectFilter(obj, predicate) {
  return Object.keys(obj)
    .filter(key => predicate(key))
    .reduce((res, key) => {
      res[key] = obj[key]
      return res
    }, {})
}

export function titleCase(str) {
  var splitStr = str.toLowerCase().split(' ')
  for (var i = 0; i < splitStr.length; i++) {
    // You do not need to check if i is larger than splitStr length, as your for does that for you
    // Assign it back to the array
    splitStr[i] = splitStr[i].charAt(0).toUpperCase() + splitStr[i].substring(1)
  }
  // Directly return the joined string
  return splitStr.join(' ')
}

export function nameList(users, expanded) {
  const names = users.map(user => titleCase(user.name)).filter(name => name !== '')
  if (names.length > 3 && !expanded) {
    return [names[0], names[names.length - 1]].join(', ') + ' et al'
  } else {
    return names.join(', ')
  }
}

export function authorList(entry, expanded) {
  if (!entry) {
    return ''
  }

  if (entry.external_db) {
    if (entry.authors?.length > 1 && expanded) {
      return `${entry.external_db} (${nameList(entry.authors)})`
    }
    return entry.external_db
  } else {
    return nameList(entry.authors || [], expanded)
  }
}

/**
 * Returns whether the used browser supports WebGL.
 * @return {boolean} Is WebGL supported.
 */
export function hasWebGLSupport() {
  let w = window
  try {
    let canvas = document.createElement('canvas')
    return !!(w.WebGLRenderingContext && (
      canvas.getContext('webgl') ||
      canvas.getContext('experimental-webgl'))
    )
  } catch (e) {
    return false
  }
}

/**
 * Returns the highest occupied energy for the given section_k_band, or null if
 * none can be found.
 *
 * For now we use section_band_gap.valence_band_max_energy as the
 * energy reference if a band gap is detected, and
 * energy_reference_fermi as the reference is a band gap is not
 * present. These represent the highest occupied energy for
 * insulators and metals respectively. Once section_results is
 * ready, it will contain the highest occupied energies explicitly.
 *
 * @param {section_k_band} section_k_band.
 * @param {scc} section_single_configuration_calculation parent section for
 * section_k_band
 *
 * @return {array} List of highest occupied energies, one for each spin
 * channel.
 */
export function getHighestOccupiedEnergy(section_k_band, scc) {
  let energyHighestOccupied
  const section_band_gap = section_k_band?.section_band_gap
  // If not energy references available, the normalization cannot be done.
  if (scc?.energy_reference_fermi === undefined && scc?.energy_reference_highest_occupied === undefined) {
    energyHighestOccupied = null
  // If a band gap is detected, it contains the most accurate highest occupied
  // energy for materials with a band gap.
  } else if (section_band_gap) {
    energyHighestOccupied = []
    for (let i = 0; i < section_band_gap.length; ++i) {
      const bg = section_band_gap[i]
      let e = bg.valence_band_max_energy
      if (e === undefined) {
        e = scc.energy_reference_fermi && scc.energy_reference_fermi[i]
        if (e === undefined) {
          e = scc.energy_reference_highest_occupied && scc.energy_reference_highest_occupied[i]
        }
      }
      energyHighestOccupied.push(e)
    }
    energyHighestOccupied = Math.max(...energyHighestOccupied)
  // Highest occupied energy reported directly by parser
  } else if (scc?.energy_reference_highest_occupied !== undefined) {
    energyHighestOccupied = scc.energy_reference_highest_occupied
  // No band gap detected -> highest occupied energy corresponds to fermi energy
  } else {
    energyHighestOccupied = Math.max(...scc.energy_reference_fermi)
  }
  return energyHighestOccupied
}

/**
 * Converts the given structure in the format used by section_results into the
 * format used by the materia-library.
 *
 * @param {object} structure.
 *
 * @return {undefined|object} If the given structure cannot be converted,
 * returns an empty object.
 */
export function toMateriaStructure(structure) {
  if (!structure) {
    return undefined
  }

  try {
    // Resolve atomic species using the labels and their mapping to chemical
    // elements.
    const speciesMap = new Map(structure.species.map(s => [s.name, s.chemical_symbols[0]]))

    const structMateria = {
      species: structure.species_at_sites.map(x => speciesMap.get(x)),
      cell: structure.lattice_vectors ? toUnitSystem(structure.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
      positions: toUnitSystem(structure.cartesian_site_positions, 'meter', {length: 'angstrom'}, false),
      fractional: false,
      pbc: structure.dimension_types ? structure.dimension_types.map((x) => !!x) : undefined
    }
    return structMateria
  } catch (error) {
    return {}
  }
}

/**
 * Given an array of numbers, calculates and returns the difference between
 * each array item and the first value in the array.
 *
 * @param {array} values Array containing values
 *
 * @return {array} Array containing the total difference values.
 */
export function diffTotal(values) {
  const diffValues = []
  const initial = values[0]
  for (let i = 0; i < values.length; i++) {
    diffValues.push(values[i] - initial)
  }
  return diffValues
}

/**
 * Formats the given number.
 *
 * @param {number} value Number to format
 * @param {decimals} decimals Number of decimals to use
 * @param {bool} scientific Whether to convert large or small values to scientific
 * form.
 *
 * @return {number} The number with new formatting
 */
export function formatNumber(value, type = 'float64', decimals = 3, scientific = true, separators = false) {
  if (isNil(value)) {
    return value
  }
  let formatted = value
  if (type?.startsWith('int')) {
    decimals = 0
  }
  if (value === 0) {
    return formatted
  }
  if (scientific) {
    if (value > 1e3 || value < 1e-3) {
      formatted = Number(Number.parseFloat(value).toExponential(decimals))
      if (separators) {
        formatted = formatted.toLocaleString()
      }
      return formatted
    }
  }
  formatted = Number(formatted.toFixed(decimals))
  if (separators) {
    formatted = formatted.toLocaleString()
  }
  return formatted
}

/**
 * Formats the given integer number.
 *
 * @param {number} value Integer to format
 *
 * @return {number} The number with new formatting
 */
export function formatInteger(value) {
  return formatNumber(value, 'int', 0, false, true)
}

/**
 * Formats the given timestamp.
 *
 * @param {number} value The timestamp to format
 * @return {str} The timestamp with new formatting
 */
export function formatTimestamp(value) {
  return value && new Date(value).toLocaleString()
}

/**
 * Determines the data type of the given metainfo.
 *
 * @param {string} quantity The metainfo name (full path). Must exist in
 * searchQuantities.json.
 *
 * @return {string} The data type of the metainfo. Can be one of the following:
 *   - number
 *   - timestamp
 *   - string
 *   - unknown
 */
export const DType = {
  Number: 'number',
  Timestamp: 'timestamp',
  String: 'string',
  Boolean: 'boolean',
  Unknown: 'unknown'
}
export function getDatatype(quantity) {
  const type = searchQuantities[quantity]?.type?.type_data

  if (isString(type) && (type.startsWith('int') || type.startsWith('float'))) {
    return DType.Number
  } else if (type === 'nomad.metainfo.metainfo._Datetime') {
    return DType.Timestamp
  } else if (type === 'str') {
    return DType.String
  } else if (type === 'bool') {
    return DType.Boolean
  } else {
    return DType.Unknown
  }
}

/**
 * Returns the unit of the given metainfo if any specified.
 *
 * @param {string} quantity The metainfo name (full path). Must exist in
 * searchQuantities.json.
 *
 * @return {string} The unit of the metainfo or undefined if no unit definition
 * found.
 */
export function getUnit(quantity) {
  return searchQuantities[quantity]?.unit
}

/**
 * Returns a function that can be used to serialize values for a given datatype.
 * @param {string} dtype The datatype
 * @param {bool} pretty Whether to serialize using a prettier, but possibly
 * lossy format.
 *
 * @return {func} Function for serializing values with the given datatype.
 */
export function getSerializer(dtype, pretty = true) {
  if (dtype === 'number') {
    return (value, units) => {
      if (isNil(value)) {
        return value
      }
      if (value instanceof Quantity) {
        if (isNil(value.value)) {
          return value.value
        }
        let label
        let valueConv
        if (units) {
          label = value.unit.label(units)
          valueConv = value.toSystem(units)
        } else {
          label = value.unit.label()
          valueConv = value.value
        }
        return `${pretty ? formatNumber(valueConv) : valueConv}${label ? ` ${label}` : ''}`
      }
      return pretty ? formatNumber(value) : value
    }
  } else if (dtype === 'timestamp') {
    return (value) => {
      if (isNil(value)) {
        return value
      }
      return pretty ? format(fromUnixTime(value / 1000), dateFormat) : value
    }
  } else {
    return (value) => value
  }
}

/**
 *
 */
export function serializeMetainfo(quantity, value, units) {
  const dtype = getDatatype(quantity)
  if (dtype === DType.Number) {
    if (!(value instanceof Quantity) && !isNil(value)) {
      const unit = getUnit(quantity) || 'dimensionless'
      value = new Quantity(value, unit)
    }
  }
  const serializer = getSerializer(dtype)
  return serializer(value, units)
}

/**
 * Returns a function that can be used to deserialize values for a given datatype.
 * @param {string} dtype The datatype
 * @param {bool} pretty Whether to deserialize using a prettier, but possibly
 * lossy format.
 *
 * @return {func} Function for deserializing values with the given datatype.
 */
export function getDeserializer(dtype, dimension) {
  if (dtype === 'number') {
    return (value, units) => {
      if (isNil(value)) {
        return value
      }
      if (isString(value)) {
        const split = value.split(' ')
        value = Number(split[0])
        if (isNaN(value)) {
          throw Error(`Could not parse the number ${split[0]}`)
        }
        return split.length === 1
          ? new Quantity(value, units?.[dimension] || 'dimensionless')
          : new Quantity(value, split[1])
      } if (isNumber(value)) {
        return new Quantity(value, units?.[dimension] || 'dimensionless')
      }
      return value
    }
  } else if (dtype === 'timestamp') {
    return (value) => {
      if (isNil(value)) {
        return value
      }
      return parseInt(value)
    }
  } else {
    return (value) => {
      const keywords = {
        true: true,
        false: false
      }
      if (value in keywords) {
        return keywords[value]
      }
      return value
    }
  }
}

/**
 * Converts a set into an array. The array will be in the insertion order of the
 * set.
 *
 * @param {Set} target Set to be converted.
 *
 * @return {number} Array containing the total difference values.
 */
export function setToArray(target) {
  if (target !== undefined && isSet(target)) {
    return [...target]
  }
  return target
}

/**
 * A simple promise based sleep.
 * @param {} ms Time in ms.
 * @returns The promise. Use .then(() => ...) to do something after sleep.
 */
export function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms))
}

/**
 * Converts a number into an approximated format with SI suffixes. The number is
 * guaranteed to be no longer than five characters, including the possible SI
 * suffix.
 *
 * @param {number} number Number to approximate
 *
 * @return {string} The number in approximated format.
 */
export function approxInteger(number) {
  // Determine the number of digits and the tier of the number
  const temp = Math.log10(Math.abs(number))
  const tier = temp / 3 | 0

  // What tier? (determines SI symbol)
  const SISymbol = ['', 'k', 'M', 'G', 'T', 'P', 'E']

  // If zero, we don't need a suffix
  if (tier === 0) return number.toFixed(0)

  // Get suffix and determine scale
  var suffix = SISymbol[tier]
  var scale = Math.pow(10, tier * 3)

  // Scale the number and count the number of decimals
  var scaled = number / scale
  let nUsed, nDecimals
  const split = scaled.toString().split('.')
  if (split.length > 1) {
    nUsed = split[0].length
    nDecimals = split[1].length
  } else {
    nUsed = split[0].length
    nDecimals = 0
  }

  // Format number and add suffix
  return scaled.toFixed(Math.min(3 - nUsed, nDecimals)) + suffix
}

/**
 * Delays the execution of the given function to the next react render cycle.
 *
 * @param {func} func Function to delay
 */
export function delay(func) {
  setTimeout(func, 0)
}

/**
 * Returns a list of linestyles.
 *
 * @param {number} nLines number of lines to plot
 */
export function getLineStyles(nLines, theme) {
  const styles = []
  const lineStyles = ['solid', 'dot', 'dashdot']
  const colors = chromaScale([theme.palette.primary.dark, theme.palette.secondary.light])
    .mode('lch').colors(nLines)
  for (let i = 0; i < nLines; ++i) {
    const line = {
      dash: lineStyles[i % lineStyles.length],
      color: colors[i],
      width: 2
    }
    styles.push(line)
  }
  return styles
}

/**
 * Returns the correct form (plural/singular) for the given word. The syntax is
 * similar to the "pluralize"-library.
 *
 * @param {string} word The word to plurarize
 * @param {count} number How many of the words exist
 * @param {boolean} inclusive Whether to prefix with the number (e.g. 3 ducks)
 * @param {boolean} format Whether to format the number.
 * @param {string} prefix An optional prefix that is placed between the number
 * and the word.
 */
export function pluralize(word, count, inclusive, format = true, prefix) {
  // Dictionary of known words. If it becomes too bothersome to keep track of
  // these words, the pluralize-library should be used instead.
  const dictionary = {
    'result': 'results',
    'entry': 'entries',
    'material': 'materials',
    'dataset': 'datasets'
  }
  const plural = dictionary[word]
  if (isNil(plural)) {
    throw Error(`The word ${word} is not in the dictionary, please add it.`)
  }
  const form = count === 1
    ? word
    : plural

  const number = inclusive
    ? format ? formatNumber(count, 'int', 0, false, true) : count
    : ''
  return `${isNil(number) ? '' : `${number} `}${isNil(prefix) ? '' : `${prefix} `}${form}`
}
