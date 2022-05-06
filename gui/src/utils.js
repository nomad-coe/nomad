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
import { cloneDeep, merge, isSet, isNil, isArray, isString, isNumber, startCase } from 'lodash'
import { Quantity } from './units'
import { format } from 'date-fns'
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
 * Map that works on n-dimensional arrays. Implemented with simple for loops for
 * performance.
 *
 * @param {*} value The original values.
 * @param {number} factor Number used for scaling.
 *
 * @return {*} A copy of the original data with numbers scaled.
 */
export function deepMap(value, func) {
  // Function for recursively applying the function
  function mapRecursive(list, newList) {
    const isScalarArray = !Array.isArray(list[0])
    if (isScalarArray) {
      for (let i = 0, size = list.length; i < size; ++i) {
        newList.push(func(list[i]))
      }
    } else {
      for (let i = 0, size = list.length; i < size; ++i) {
        const iList = []
        newList.push(iList)
        mapRecursive(list[i], iList)
      }
    }
  }
  // If given a scalar variable, simply apply the function. If a list is given,
  // perform the mapping recursively.
  let newValue
  if (!isArray(value)) {
    newValue = func(value)
  } else {
    newValue = []
    mapRecursive(value, newValue)
  }
  return newValue
}

/**
 * Scales values of an n-dimensional array.
 *
 * @param {n-dimensional array} value The original values.
 * @param {number} factor Scaling factor.
 *
 * @return {n-dimensional array} A copy of the original data with numbers scaled.
 */
export function scale(value, factor) {
  return deepMap(value, x => x * factor)
}

/**
 * Adds the given value to all elements of an n-dimensional array.
 *
 * @param {n-dimensional array} value The original values.
 * @param {number} addition Value to add.
 *
 * @return {n-dimensional array} A copy of the original data with the addition.
 */
export function add(value, addition) {
  return deepMap(value, x => x + addition)
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
      cell: structure.lattice_vectors
        ? new Quantity(structure.lattice_vectors, 'meter').to('angstrom').value()
        : undefined,
      positions: new Quantity(structure.cartesian_site_positions, 'meter').to('angstrom').value(),
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
 * @param {string} type Number data type.
 * @param {number} decimals Number of decimals to use. Note
 * @param {string} mode The formatting mode to use. One of: 'scientific', 'separators', 'standard'
 * available for scientific formatting.
 *
 * @return {string} The number as a string with new formatting
 */
export function formatNumber(value, type = DType.Float, mode = 'scientific', decimals = 3) {
  // Nill values
  if (isNil(value)) {
    return value
  }
  // Zero
  if (value === 0) {
    return value.toString()
  }

  // Scientific format. parseFloat is used to get rid of trailig insignificant
  // zeros.
  if (mode === 'scientific') {
    const absValue = Math.abs(value)
    return (absValue > 1e3 || absValue < 1e-2)
      ? parseFloat(value.toExponential(decimals)).toExponential()
      : parseFloat(Number(value.toFixed(decimals))).toString()
  } else {
    // Integers in-non-scientific formats are always shown without decimals
    if (type?.startsWith('int')) {
      decimals = 0
    }
    const formatted = Number(value.toFixed(decimals))
    // Format with separators
    if (mode === 'separators') {
      return formatted.toLocaleString()
    // Standard formatting
    } else if (mode === 'standard') {
      return formatted.toString()
    // SI formatting
    } else if (mode === 'SI') {
      return approxInteger(formatted)
    } else {
      throw Error('Unknown formatting mode')
    }
  }
}

/**
 * Formats the given integer number.
 *
 * @param {number} value Integer to format
 *
 * @return {number} The number with new formatting
 */
export function formatInteger(value) {
  return formatNumber(value, DType.Int, 'separators', 0)
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
  Int: 'int',
  Float: 'float',
  Timestamp: 'timestamp',
  String: 'string',
  Enum: 'enum',
  Boolean: 'boolean',
  Unknown: 'unknown'
}
const numericTypes = new Set([DType.Int, DType.Float])
export function getDatatype(quantity) {
  const type_data = searchQuantities[quantity]?.type?.type_data
  const type_kind = searchQuantities[quantity]?.type?.type_kind

  if (isString(type_data) && type_data.startsWith('int')) {
    return DType.Int
  } else if (isString(type_data) && type_data.startsWith('float')) {
    return DType.Float
  } else if (type_data === 'nomad.metainfo.metainfo._Datetime') {
    return DType.Timestamp
  } else if (type_data === 'str') {
    return DType.String
  } else if (type_kind === 'Enum') {
    return DType.Enum
  } else if (type_data === 'bool') {
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
  if (numericTypes.has(dtype)) {
    return (value, units) => {
      if (isNil(value)) {
        return value
      }
      if (value instanceof Quantity) {
        const newVal = units ? value.toSystem(units) : value
        const label = newVal.label()
        const valueConv = newVal.value()
        return `${pretty ? formatNumber(valueConv) : valueConv}${label ? ` ${label}` : ''}`
      }
      return pretty ? formatNumber(value) : value
    }
  } else if (dtype === DType.Timestamp) {
    return (value) => {
      if (isNil(value)) {
        return value
      }
      return pretty ? format(new Date(value), dateFormat) : value
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
  if (dtype === DType.Int || dtype === DType.Float) {
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
  if (numericTypes.has(dtype)) {
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
  } else if (dtype === DType.Timestamp) {
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
    'search result': 'search results',
    'entry': 'entries',
    'value': 'values',
    'material': 'materials',
    'dataset': 'datasets',
    'item': 'items'
  }
  const plural = dictionary[word]
  if (isNil(plural)) {
    throw Error(`The word ${word} is not in the dictionary, please add it.`)
  }
  const form = count === 1
    ? word
    : plural

  const number = inclusive
    ? format ? formatNumber(count, DType.Int, 'separators', 0) : count
    : ''
  return `${isNil(number)
    ? ''
    : `${number} `}${isNil(prefix)
    ? ''
    : `${prefix} `}${form}`
}

/**
 * Used to create a label from the metainfo name.
 * @param {str} name Metainfo name
 * @returns A formatted label constructed from the metainfo name.
 */
export function getLabel(name) {
  let label = searchQuantities[name]?.name || name
  label = label.replace(/_/g, ' ')
  label = startCase(label)
  return label
}

/**
 * Used for testing purposes: setting a data-testid to this value signals that the
 * component waits for further updates of some kind. This is used by some automated tests
 * to determine if a component is *really* done rendering (when we don't really care what
 * we are waiting for, and just want a way to determine when we're done). When done,
 * the data-testid should be removed or updated to something else.
 */
export const isWaitingForUpdateTestId = 'waiting-for-update'
