import { unit } from 'mathjs'
import { cloneDeep, merge } from 'lodash'

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
 * Used to convert numeric values from one unit to another. Works on
 * n-dimensional arrays and implemented as a relatively simple for loop for
 * performance. If conversion times become an issue, it might be worthwhile to
 * look at vectorization with WebAssembly.
 *
 * @param {*} value The values to convert
 * @param {*} from Original unit.
 * @param {*} to Target unit.
 *
 * @return {*} A copy of the original data with units converted.
 */
export function convert(value, from, to) {
  // Determine the scaling factor
  let factor = unit(1, from).toNumber(to)

  // Convert arrays
  function scaleRecursive(list, newList) {
    let isScalarArray = !Array.isArray(list[0])
    if (isScalarArray) {
      for (let i = 0, size = list.length; i < size; ++i) {
        newList.push(list[i] * factor)
      }
    } else {
      for (let i = 0, size = list.length; i < size; ++i) {
        let iList = []
        newList.push(iList)
        scaleRecursive(list[i], iList)
      }
    }
  }
  let isArray = Array.isArray(value)
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
