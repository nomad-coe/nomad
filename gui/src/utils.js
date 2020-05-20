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
