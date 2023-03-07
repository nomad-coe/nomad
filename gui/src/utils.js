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
import { cloneDeep, merge, isSet, isNil, isArray, isString, isNumber, isPlainObject, startCase } from 'lodash'
import { Quantity } from './units'
import { format } from 'date-fns'
import { dateFormat, guiBase, apiBase } from './config'
import searchQuantities from './searchQuantities.json'
import parserMetadata from './parserMetadata.json'
const crypto = require('crypto')

const parserLabels = Object.keys(parserMetadata).reduce((result, key) => {
  const parser = parserMetadata[key]
  result[`parsers/${parser.codeName}`] = parser.codeLabel
  result[parser.codeName] = parser.codeLabel
  return result
}, {})

export const isEquivalent = (a, b) => {
  // Create arrays of property names
  const aProps = Object.getOwnPropertyNames(a)
  const bProps = Object.getOwnPropertyNames(b)

  // If number of properties is different,
  // objects are not equivalent
  if (aProps.length !== bProps.length) {
    return false
  }

  for (let i = 0; i < aProps.length; i++) {
    const propName = aProps[i]

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
 * Recursive object traversal.
 *
 * @param {*} data The data to traverse.
 * @param {bool} fullPath Whether to return the full path as a list of keys. Defaults to
 *   false which means that only the last key is returned.
 * @return List of pairs of keys and values.
 */
export function traverseDeep(object, fullPath = false) {
  const items = []
  const traverse = (obj, path) => {
    for (const [key, value] of Object.entries(obj)) {
      const newPath = [...path]
      newPath.push(key)
      if (value && isPlainObject(value)) {
        traverse(value, newPath)
      } else {
        items.push([fullPath ? [...newPath] : key, value])
      }
    }
  }
  traverse(object, [])
  return items
}

/**
 * Used to resolve a property from a nested object. Also supports indices in the
 * path.
 *
 * @param {object} data The data to traverse.
 * @param {str} path Path to use
 * @return The value at the given path or undefined if none can be found.
 */
export function getDeep(data, path, separator = '.') {
    const segments = path.split(separator)
    return segments.reduce((current, segment) => current && current[segment], data)
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
export function mapDeep(value, func) {
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
  return mapDeep(value, x => x * factor)
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
  return mapDeep(value, x => x + addition)
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
  const cloned = cloneDeep(target)
  const val = merge(cloned, source)
  return val
}

export function arraysEqual(_arr1, _arr2) {
  if (!Array.isArray(_arr1) || !Array.isArray(_arr2) || _arr1.length !== _arr2.length) {
    return false
  }

  const arr1 = _arr1.concat().sort()
  const arr2 = _arr2.concat().sort()

  for (let i = 0; i < arr1.length; i++) {
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
  const splitStr = str.toLowerCase().split(/[_ ]/)
  for (let i = 0; i < splitStr.length; i++) {
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
 * @param {*} entry metadata retrived from NOMAD api
 * @returns The name from the metadata, a fallback entry name for simulations, or null
 */
export function entryName(entry) {
  const results = entry?.results
  const name = entry?.entry_name
  if (name) {
    return name
  }
  const type = entryType(entry)
  if (!type) {
    return null
  }
  const formula = results?.material?.chemical_formula_hill
  if (formula) {
    return `${formula} ${type}`
  }
  return type
}

/**
 * @param {*} entry metadata retrived from NOMAD api
 * @returns The type from the metadata, a fallback entry type for simulations, or null
 */
export function entryType(entry) {
  const results = entry?.results
  let type = entry?.entry_type
  if (type) {
    return type
  }
  if (!(results?.method?.simulation && results?.material)) {
    const parser = entry?.parser_name
    if (!parser) {
      return null
    }
    if (parser === 'parsers/eels') {
      return 'EELS DB entry'
    }
    const program = parserLabels[parser]
    if (!program) {
      return null
    }
    return `${program} simulation`
  }
  type = ''
  let workflow = results.method.workflow_name?.[0]
  if (workflow === 'single_point') {
    workflow = null
  }
  if (workflow) {
    type += workflow + ' '
  }
  const program = results.method.simulation.program_name
  if (program) {
    type += program + ' '
  }
  if (type === '') {
    return null
  }
  type += 'simulation'
  return type
}

/**
 * Returns whether the used browser supports WebGL.
 * @return {boolean} Is WebGL supported.
 */
let hasWebGL
export function hasWebGLSupport() {
  if (isNil(hasWebGL)) {
    const w = window
    try {
      const canvas = document.createElement('canvas')
      hasWebGL = !!(w.WebGLRenderingContext && (
        canvas.getContext('webgl') ||
        canvas.getContext('experimental-webgl'))
      )
    } catch (e) {
      hasWebGL = false
    }
  }
  return hasWebGL
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
 * @param {string} mode The formatting mode to use. One of:
 *  - 'scientific'
 *  - 'separators'
 *  - 'standard'
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
  if (isNil(value)) {
    return value
  }
  if (value.search(/([+-][0-9]{2}:[0-9]{2}|Z)\b/) === -1) { // search for timezone information
    try {
      // assume UTC timestamp from server and attempt to manually add UTC timezone,
      // new Date will wrongly assume local timezone.
      const result = new Date(`${value}Z`).toLocaleString()
      if (result !== 'Invalid Date') {
        return result
      }
    } catch {}
  }
  // create string for timestamp with correct timezone information or fallback
  return new Date(value).toLocaleString()
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
  if (typeof quantity === 'string') {
    quantity = searchQuantities[quantity]
  }
  const type_data = quantity?.type?.type_data
  const type_kind = quantity?.type?.type_kind

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
  const suffix = SISymbol[tier]
  const scale = Math.pow(10, tier * 3)

  // Scale the number and count the number of decimals
  const scaled = number / scale
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
 * Forces the given function to be placed in the Javascript event queue and
 * executed once the queue is free. Returns a promise that also correctly
 * catches errors within the execution.
 *
 * @param {func} func Function to delay
 */
export function delay(func) {
  return new Promise((resolve, reject) => {
    setTimeout(() => {
      try {
        func()
        resolve()
      } catch (e) {
        reject(e)
      }
    })
  })
}

/**
 * Returns the correct form (plural/singular) for the given word. The syntax is
 * similar to the "pluralize"-library.
 *
 * @param {string} word The word to pluralize
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
    'item': 'items',
    'upload': 'uploads',
    'code': 'codes',
    'manager': 'managers',
    'filter': 'filters'
  }
  const words = word.trim().split(" ")
  const lastWord = words[words.length - 1]
  let plural = dictionary[lastWord]
  if (isNil(plural)) {
    throw Error(`The word ${word} is not in the dictionary, please add it.`)
  }
  words[words.length - 1] = plural
  plural = words.join(" ")
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
 * Used to create a formatted label for a metainfo name or value. Replaces
 * underscores with whitespace and capitalizes the first letters.
 *
 * @param {str} name Metainfo name
 * @returns A formatted label constructed from the metainfo name.
 */
export function formatLabel(label) {
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

/**
 * Removes indentation from the given multiline code string
 * @param {*} code
 */
export function stripIndent(code) {
  let lines = code.split('\n')
  const start = lines.findIndex(line => line.trim() !== '')
  lines = lines.slice(start)
  if (lines.length === 0) {
    return ''
  }
  const indent = lines[0].length - lines[0].trimStart().length
  return lines.map(line => line.substring(indent)).join('\n')
}

/**
 * Used to imperatively get a path that works for react-router and points to the
 * currently displayed page. The useLocation hook should be preferred, but
 * sometimes getting the path imperatively can prevent unnecessary renders.
 *
 * @returns A formatted label constructed from the metainfo name.
 */
export function getLocation() {
  return `${window.location.pathname.slice(guiBase.length)}`
}

/**
 * Calculates the entryId given the uploadId, the mainfile, and, optionally, the mainfileKey.
 */
 export function generateEntryId(uploadId, mainfile, mainfileKey) {
  const hash = crypto.createHash('sha512').update(uploadId).update(mainfile)
  if (mainfileKey) {
    hash.update(mainfileKey)
  }
  return hash.digest('base64').slice(0, 28)
    .replace(/\+/g, '-') // Convert '+' to '-'
    .replace(/\//g, '_') // Convert '/' to '_'
}

/**
 * Utilities parse and analyze *nomad urls*. A nomad url identifies some *base resource*,
 * which can be:
 *  1) a nomad deployment,
 *  2) an upload,
 *  3) an archive, or
 *  4) a metainfo schema.
 *
 * A *metainfo schema* is a json data structure, which is either
 *  A) the system metainfo, generated from the python classes,
 *  B) provided in an entry uploaded by a user, or
 *  C) a *frozen* metainfo schema.
 * A frozen metainfo schema is a historical "snapshot" of a metainfo schema provided using
 * method A or B. The only difference between the frozen version and the original version is that
 * in the frozen version, references to other metainfo schemas are replaced with references to
 * frozen versions of these schemas. The purpose of this is to allow us to handle changing
 * metainfo schemas.
 *
 * In addition to specifying the base resource, a nomad url can specify a location *within*
 * this resource (identifying a raw file or folder in an upload, a location in an archive, or
 * a section in a metainfo schema). For urls identifying sections in a metainfo schema,
 * it is generally also possible to specify a versionHash, which uniquely identifies a frozen
 * version of this section definition.
 *
 * Nomad urls can be absolute or relative. Absolute urls contain all information needed to
 * locate the resource, including the nomad deployment url. Relative urls can only be resolved
 * when given a *base url*, which defines our "point of origin".
 *
 * Possible formats (expressions in [brackets] are optional):
 * ----------------------------------------------------------
 *  Absolute urls:
 *    <deploymentUrl>
 *    <deploymentUrl>/uploads/<uploadId> [ /raw/<rawPath> [ #<dataPath> [ @<versionHash> ] ] ]
 *    <deploymentUrl>/uploads/<uploadId>/archive/<entryid> [ #<dataPath> [ @<versionHash> ] ]
 *    <deploymentUrl>/uploads/<uploadId>/archive/mainfile/<mainfile> [ #<dataPath> [ @<versionHash> ] ]
 *  Urls relative to the current deployment:
 *    ../uploads/<uploadId> [ /raw/<rawPath> [ #<dataPath> [ @<versionHash> ] ] ]
 *    ../uploads/<uploadId>/archive/<entryid> [ #<dataPath> [ @<versionHash> ] ]
 *    ../uploads/<uploadId>/archive/mainfile/<mainfile> [ #<dataPath> [ @<versionHash> ] ]
 *    <qualifiedName> (TODO: how to handle versions, deploymentUrl etc)
 *  Urls relative to the current upload:
 *    ../upload [ /raw/<rawPath> [ #<dataPath> [ @<versionHash> ] ] ]
 *    ../upload/archive/<entryid> [ #<dataPath> [ @<versionHash> ] ]
 *    ../upload/archive/mainfile/<mainfile> [ #<dataPath> [ @<versionHash> ] ]
 *  Urls relative to the current data (i.e. current archive or metainfo schema json):
 *    #<dataPath> (preferred)
 *    /<dataPath>
 * Note:
 *  - An <deploymentUrl> is a normal url, starting with "http://" or "https://", locating
 *    the api of the nomad deployment. Should always end with "/api".
 *    Example:
 *      https://nomad-lab.eu/prod/rae/api
 *  - Urls with versionHash are considered to be relative to the deployment rather than the
 *    upload or data.
 *  - The rawPath and mainFile paths need to be escaped with urlEncodePath to ensure a valid url.
 *    (i.e. each segment needs to individually be escaped using encodeURIComponent)
 *  - The dataPath, if provided, must always start from the root of the archive.
 *    Leading and duplicate slashes are ignored.
 *  - If the first segment in dataPath is 'definitions' or 'packages', the url is recognized
 *    as a reference to a metainfo schema.
 *  - Metainfo sections defined in python can be referred using the qualifiedName, which
 *    is a "python style" name of alphanumerics separated by dots.
 *  - If no versionHash is specified for a url identifying a metainfo schema, it means
 *    we refer to the version defined by the base url (if the url is schema-relative), otherwise
 *    the latest version (in the nomad deployment in question).
 */

export const systemMetainfoUrl = `system-metainfo` // TODO: should use absolute url with hash when implemented

/**
 * Enum for the `type` attribute of parsed/resolved url objects
 */
export const refType = Object.freeze({
  deployment: 'deployment',
  upload: 'upload',
  archive: 'archive',
  metainfo: 'metainfo'
})

/**
 * Enum for the `relativeTo` attribute of parsed/resolved url objects
 */
export const refRelativeTo = Object.freeze({
  deployment: 'deployment',
  upload: 'upload',
  data: 'data',
  nothing: null
})

/**
 * Parses a nomad url and returns an object with the following attributes:
 *  url
 *    the original url string.
 *  type
 *    One of: refType.deployment | refType.upload | refType.archive | refType.metainfo
 *  relativeTo
 *    One of: refRelativeTo.deployment | refRelativeTo.upload | refRelativeTo.data | refRelativeTo.nothing (null)
 *  deploymentUrl
 *    The nomad deployment url, if it can be determined.
 *  uploadId
 *    The uploadId, if it can be determined.
 *  entryId
 *    The entryId, if it can be determined.
 *  mainfile
 *    The mainfile path, if it can be determined (unescaped!).
 *  path
 *    The path within the base resource (dataPath | rawPath (unescaped!)), if specified.
 *    If it is a dataPath, it will have a leading slash.
 *  qualifiedName
 *    The qualifiedName, if specified.
 *  versionHash
 *    The versionHash, if specified.
 *  isResolved
 *    True if the url is fully resolved, i.e. we have everything we need to fetch the data.
 *  isExternal
 *    If the url refers to a resource in an external nomad deployment. Note, for relative
 *    urls the value will be undefined, which means we don't know.
 *
 * If the url cannot be parsed, an error is thrown.
 *
 * @param {*} url either a string, or an object. If an object, it is expected to be parsed already,
 *    and we will return the object "as is" directly, without any parsing.
 */
export function parseNomadUrl(url) {
  if (isPlainObject(url) && url.url && url.type) return url // already parsed, pass through.
  const prefix = `Could not parse nomad url "${url}": `
  if (!url) throw new Error(prefix + 'empty value')
  if (typeof url !== 'string') throw new Error(prefix + 'bad type, expected string, got ' + typeof url)

  if (url === systemMetainfoUrl) {
    // TODO proper handling, using a url containing deploymentUrl + versionHash somehow?
    return {
      url,
      relativeTo: refRelativeTo.deployment,
      type: refType.metainfo,
      deploymentUrl: apiBase,
      isResolved: true,
      isExternal: false
    }
  }

  let relativeTo, type, deploymentUrl, uploadId, entryId, mainfile, path, qualifiedName, versionHash
  let dataPath, rest, rawPath
  if (url.startsWith('/')) {
    dataPath = url
    rest = ''
  } else {
    const bracketPos = url.indexOf('#')
    dataPath = bracketPos !== -1 ? url.slice(bracketPos + 1) : undefined
    rest = bracketPos !== -1 ? url.slice(0, bracketPos) : url
  }

  if (rest.startsWith('http://') || rest.startsWith('https://')) {
    // Url includes deploymentUrl
    let apiPos = rest.indexOf('/api/')
    if (apiPos === -1 && rest.endsWith('/api')) apiPos = rest.length - 4
    if (apiPos === -1) throw new Error(prefix + 'absolute nomad deployment url does not contain "/api"')
    deploymentUrl = url.slice(0, apiPos + 4)
    rest = rest.slice(apiPos + 5)
    if (rest && !rest.startsWith('uploads/')) throw new Error(prefix + 'expected "/uploads/<uploadId>" in absolute url')
    relativeTo = null
  } else if (url.startsWith('../')) {
    rest = rest.slice(3)
  } else if (url.startsWith('#') || url.startsWith('/')) {
    relativeTo = refRelativeTo.data
  } else if (url.match(/^([a-zA-Z]\w*\.)*[a-zA-Z]\w*$/)) {
    qualifiedName = url
    relativeTo = refRelativeTo.deployment
  } else {
    throw new Error(prefix + 'bad url')
  }
  const restParts = rest.split('/')
  if ((deploymentUrl && rest) || url.startsWith('../')) {
    // Expect upload ref
    if (restParts[0] === 'uploads') {
      // Ref with upload id
      if (!deploymentUrl) {
        relativeTo = refRelativeTo.deployment
      }
      if (restParts.length === 1) throw new Error(prefix + 'expected "/uploads/<uploadId>" in url')
      uploadId = restParts[1]
      restParts.splice(0, 2)
    } else if (restParts[0] === 'upload') {
      // Relative ref
      relativeTo = refRelativeTo.upload
      restParts.splice(0, 1)
    } else {
      throw new Error(prefix + 'expected "/upload" or "/uploads" in url')
    }
    if (restParts.length) {
      // There is more. Expect "raw" or "archive"
      if (restParts[0] === 'raw') {
        rawPath = restParts.slice(1).map(decodeURIComponent).join('/')
      } else if (restParts[0] === 'archive') {
        if (restParts.length === 1) throw new Error(prefix + '"archive" must be followed by entry id or "mainfile"')
        if (restParts[1] === 'mainfile') {
          if (restParts.length === 2) throw new Error(prefix + '"mainfile" must be followed by a mainfile path')
          mainfile = restParts.slice(2).map(decodeURIComponent).join('/')
        } else {
          if (restParts.length !== 2) throw new Error(prefix + 'unexpected path element after entry id')
          entryId = restParts[1]
        }
      } else {
        throw new Error(prefix + 'expected "raw" or "archive" after upload ref')
      }
    }
  }
  if (qualifiedName) {
    // Refers to the global schema
    type = refType.metainfo
  } else if (deploymentUrl && !uploadId) {
    // Pure deployment url
    type = refType.deployment
  } else if (dataPath !== undefined) {
    // Has dataPath
    if (!url.startsWith('#') && !url.startsWith('/') && !entryId && !mainfile && !rawPath) throw new Error(prefix + 'Unexpected dataPath without entry reference')
    mainfile = mainfile || rawPath
    const atPos = dataPath.indexOf('@')
    if (atPos !== -1) {
      versionHash = dataPath.slice(atPos + 1)
      dataPath = dataPath.slice(0, atPos)
    }
    const dataPathSegments = dataPath.split('/').filter(segment => segment)
    const firstDataPathSegment = dataPathSegments[0]
    type = (firstDataPathSegment === 'definitions' || firstDataPathSegment === 'packages') ? refType.metainfo : refType.archive
    dataPath = '/' + dataPathSegments.join('/') // Ensure leading slash, but no duplicate slashes
    path = dataPath
  } else if (entryId || mainfile) {
    // Refers to an archive, but has no dataPath
    type = refType.archive
  } else {
    // Refers to an upload
    type = refType.upload
    path = rawPath
  }
  if (versionHash !== undefined) {
    if (type !== refType.metainfo) throw new Error(prefix + 'versionHash can only be specified for metainfo urls.')
    if (!versionHash.match(/\w+/)) throw new Error(prefix + 'bad versionHash provided')
    if (relativeTo) {
      relativeTo = refRelativeTo.deployment
    }
  }

  if (uploadId && mainfile) {
    entryId = generateEntryId(uploadId, mainfile)
  }

  return {
    url,
    relativeTo,
    type,
    deploymentUrl,
    uploadId,
    entryId,
    mainfile,
    path,
    qualifiedName,
    versionHash,
    isResolved: !relativeTo,
    isExternal: deploymentUrl ? (deploymentUrl !== apiBase) : undefined,
    toString: () => { return url + ' (parsed)' }
  }
}

/**
 * Resolves the url with respect to the given baseUrl. Each of the two url arguments can
 * be either a string or a parsed/resolved nomad url object (i.e. an object returned by calling
 * parseNomadUrl or resolveNomadUrl). If the url is absolute, no baseUrl is required and the argument
 * is ignored. If url is relative, however, baseUrl is required, and it needs to be either a string
 * which is an absolute url, or a resolved nomad url object.
 *
 * If succesful, returns an resolved nomad url object, which is similar to the one returned
 * by parseNomadUrl, but
 *  1) It additionally has the attribute baseUrl, which stores the baseUrl argument,
 *  2) the isResolved flag should be guaranteed to be true, and
 *  3) all applicable attributes (like deploymentUrl, uploadId, entryId, etc) are set
 *     to the resolved valules.
 *
 * NOTE, if you pass an object as the url argument, and it is not already resolved, the
 * original object will be modified by this call, and the modified version is also the
 * object that will be returned by the function.
 */
export function resolveNomadUrl(url, baseUrl) {
  const parsedUrl = parseNomadUrl(url)
  if (parsedUrl.isResolved) return parsedUrl
  parsedUrl.baseUrl = baseUrl
  const prefix = `Could not normalize url "${parsedUrl.url}" with baseUrl "${baseUrl?.url || baseUrl}": `

  if (parsedUrl.relativeTo) {
    // Url is relative.
    if (!baseUrl) throw new Error(prefix + 'a baseUrl is required.')
    const parsedBaseUrl = parseNomadUrl(baseUrl)
    if (!parsedBaseUrl.isResolved) throw new Error(prefix + 'unresolved baseUrl')
    // Copy data from parsedBaseUrl
    parsedUrl.deploymentUrl = parsedBaseUrl.deploymentUrl // Should always be copied
    if (parsedUrl.relativeTo === refRelativeTo.upload) {
      if (!parsedBaseUrl.uploadId) throw new Error(prefix + 'missing information about uploadId')
      parsedUrl.uploadId = parsedBaseUrl.uploadId
      if (parsedUrl.mainfile) {
        parsedUrl.entryId = generateEntryId(parsedUrl.uploadId, parsedUrl.mainfile)
      }
    } else if (parsedUrl.relativeTo === refRelativeTo.data) {
      if (!parsedBaseUrl.entryId || !parsedBaseUrl.uploadId) throw new Error(prefix + 'missing information about entryId')
      parsedUrl.uploadId = parsedBaseUrl.uploadId
      parsedUrl.entryId = parsedBaseUrl.entryId
    }
  }

  parsedUrl.isExternal = !!parsedUrl.deploymentUrl && parsedUrl.deploymentUrl !== apiBase
  parsedUrl.isResolved = true
  return parsedUrl
}

/**
 * Works like resolveNomadUrl, but does not throw errors in case of failure. Instead:
 *  1) If url is empty, we return null
 *  2) If the url is non-empty and unresolvable, we return an *resolve failed object*, which
 *     is an object containing {url, baseUrl, error}
 */
export function resolveNomadUrlNoThrow(url, baseUrl) {
  if (!url) return null
  try {
    return resolveNomadUrl(url, baseUrl)
  } catch (error) {
    return {url, baseUrl, error}
  }
}

/**
 * Relativizes a url with respect to the provided deploymentUrl, uploadId, entryId, and
 * returns the shortest possible relative url, as a string.
 * This method thus basically does the *opposite* of resolveNomadUrl. You can specify either none,
 * the first, the two first, or all three arguments of the (deploymentUrl, uploadId, entryId)
 * tuple. The url provided must be absolute.
 */
export function relativizeNomadUrl(url, deploymentUrl, uploadId, entryId) {
  const parsedUrl = parseNomadUrl(url)
  if (!parsedUrl.isResolved) {
    throw new Error(`Absolute url required, got ${url}.`)
  }
  if (!parsedUrl.uploadId) {
    throw new Error(`Expected url to specify an upload, got ${url}`)
  }
  if (parsedUrl.deploymentUrl !== deploymentUrl) {
    // Nothing to relativize
    return normalizeNomadUrl(parsedUrl)
  }
  if (parsedUrl.entryId === entryId) {
    // Same deployment and entry
    return '#' + (parsedUrl.path || '/')
  }
  // Same deployment, possibly also same upload
  let rv = parsedUrl.uploadId === uploadId ? '../upload' : `../uploads/${parsedUrl.uploadId}`
  if (parsedUrl.type === refType.archive || parsedUrl.type === refType.metainfo) {
    rv = `${rv}/archive/${parsedUrl.entryId}`
    if (parsedUrl.path) {
      rv = `${rv}#${parsedUrl.path}`
    }
  } else if (parsedUrl.type === refType.upload) {
    if (parsedUrl.path) {
      rv = `${rv}/raw/${parsedUrl.path}`
    }
  }
  return rv
}

/**
 * Returns a url string which is "normalized", i.e. it is absolute and have the preferred
 * url form for referencing the resource/data. Normalized urls always prefer identifying
 * entries with entryId rather than by specifying a mainfile.
 * @param {*} url A string or resolved nomad url object.
 */
export function normalizeNomadUrl(url) {
  const parsedUrl = parseNomadUrl(url)
  if (!parsedUrl.isResolved) throw new Error(`Failed to normalize url ${url.url}: provided url is unresolved`)
  if (parsedUrl.type === refType.deployment) {
    return parsedUrl.deploymentUrl
  } else if (parsedUrl.type === refType.upload) {
    return `${parsedUrl.deploymentUrl}/uploads/${parsedUrl.uploadId}` + (parsedUrl.path ? '/raw/' + parsedUrl.path : '')
  } else if (parsedUrl.entryId) {
    return `${parsedUrl.deploymentUrl}/uploads/${parsedUrl.uploadId}/archive/${parsedUrl.entryId}` + (
      parsedUrl.path ? '#' + parsedUrl.path + (parsedUrl.versionHash ? '@' + parsedUrl.versionHash : '') : '')
  } else if (parsedUrl.qualifiedName) {
    return parsedUrl.qualifiedName
  }
  throw new Error(`Failed to normalize url ${url.url}: unhandled url format`) // Should not happen
}

/**
 * Utility for creating an absolute upload url, given deploymentUrl, uploadId and an UNESCAPED rawPath
 */
export function createUploadUrl(deploymentUrl, uploadId, rawPathUnescaped) {
  const rawPathEscaped = urlEncodePath(rawPathUnescaped || '')
  return `${deploymentUrl}/uploads/${uploadId}/raw/${rawPathEscaped}`
}

/**
 * Utility for creating an absolute entry url, given deploymentUrl, uploadId, entryId, and dataPath
 */
export function createEntryUrl(deploymentUrl, uploadId, entryId, dataPath) {
  let rv = `${deploymentUrl}/uploads/${uploadId}/archive/${entryId}`
  if (dataPath) {
    rv += '#' + (!dataPath.startsWith('/') ? '/' : '') + dataPath
  }
  return rv
}

/**
 * Given parsedUrl, which should be a resolved, parsed url, specifying an entry and some dataPath
 * within this entry (the path may be empty), append the dataPathSuffix, and return the resulting
 * url as a string. If the dataPath of parsedUrl is empty, it will be interpreted to mean the archive
 * root element. The dataPathSuffix should not start with a '/'.
 */
export function appendDataUrl(parsedUrl, dataPathSuffix) {
  if (!parsedUrl.isResolved) throw new Error(`appendDataUrl: a resolved url required, got ${parsedUrl}`)
  if (parsedUrl.type !== refType.archive) throw new Error(`appendDataUrl: an archive url required, got ${parsedUrl}`)
  return urlJoin(
    `${parsedUrl.deploymentUrl}/uploads/${parsedUrl.uploadId}/archive/${parsedUrl.entryId}#${parsedUrl.path || '/'}`,
    dataPathSuffix)
}

/**
 * A simple method for joining url components. The componenst should be strings, and the return
 * value is also a string. A '/' is inserted between all components, if needed (i.e. if the preceding
 * component does not end with a '/').
 *
 * Note that the standard join method defined in the path library does not quite work like
 * this, it reduces all occurances of '//' to '/' in the returned value, and we don't want this
 * for urls (we want to keep the double slash in expressions like 'http://' etc).
 */
export function urlJoin(...components) {
  let rv = ''
  for (const component of components) {
    if (component) {
      if (rv && !rv.endsWith('/')) {
        rv += '/'
      }
      rv += component
    }
  }
  return rv
}

/**
 * Encodes a file system path (which may include various characters not allowed in urls, or that
 * have special meanings in urls, like ' ' and '#' and '?') into a path that is url safe.
 */
export function urlEncodePath(path) {
  return path.split('/').map(encodeURIComponent).join('/')
}

/**
 * Decodes a url encoded file system path (the inverse of urlEncodePath)
 */
export function urlDecodePath(urlPath) {
  return urlPath.split('/').map(decodeURIComponent).join('/')
}

/**
 * Given an internal url and a data object, gets the value at the corresponding location (path)
 * in this data object. The data object should be an archive or a metainfo data object. The
 * url parameter can be a single url or an array of urls. If it is an array of urls, the method
 * will return an array of values instead of a single value. The url argument can be empty,
 * in which case null will be returned. Otherwise, the url(s) MUST be internal (i.e. data-relative,
 * starting with a '#' or a '/'), and they can be either strings or parsed nomad url objects.
 */
export function resolveInternalRef(url, data, throwOnError = true) {
  if (!url) return null
  if (Array.isArray(url)) {
    return url.map(ref => resolveInternalRefSingle(ref, data, throwOnError))
  }
  return resolveInternalRefSingle(url, data, throwOnError)
}

function resolveInternalRefSingle(url, data) {
  let path
  if (url.relativeTo === refRelativeTo.data) {
    // We have a parsed url obj
    path = url.path
  } else {
    // Should be a string object
    if (typeof url !== 'string' || (!url.startsWith('#') && !url.startsWith('/'))) throw new Error(`Expected internal ref, got ${url}`)
    path = url.slice(1)
  }
  const segments = path.split('/').filter(segment => segment !== '')
  let current = data
  for (const segment of segments) {
    if (!current) throw new Error(`Path does not exist: ${path}`)
    if (isNaN(segment)) {
      current = current[segment] || current._properties?.[segment] || current.sub_section?._properties?.[segment]
    } else {
      if (current.repeats && current.m_def) {
      } else {
        let index = parseInt(segment)
        if (index < 0) index = current.length - index
          current = current[index]
      }
    }
  }
  return current
}

/**
 * Defines an alphabet in capital letters
 */
export const alphabet = [
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
  'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'
]
