import { unit } from 'mathjs'

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
 * look at vectorization with SIMD.
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
