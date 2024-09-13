// Returns groups in arbitrary order without duplicates with the given IDs
export async function fetchGroupsByIds(api, groupIds) {
  if (groupIds.length === 0) {
    return []
  }

  const uniqueIds = [...new Set(groupIds)]
  return new Promise((resolve, reject) => {
    api.getGroups({group_id: uniqueIds})
      .then(resolve)
      .catch(error => reject(new Error(`Unable to fetch the groups by IDs (${error})`)))
  })
}

// Returns users in arbitrary order without duplicates with the given IDs
export async function fetchUsersByIds(api, userIds) {
  if (userIds.length === 0) {
    return []
  }

  const uniqueIds = [...new Set(userIds)]
  return new Promise((resolve, reject) => {
    api.getUsers({user_id: uniqueIds})
      .then(resolve)
      .catch(error => reject(new Error(`Unable to fetch the users by IDs (${error})`)))
  })
}

/**
  * Returns groups in arbitrary order without duplicates whose names contain the search terms
  * {string} search_terms - space-separated substrings, all must be present in the group name (no case)
  */
export async function fetchGroupsByNameSearch(api, search_terms) {
  if (search_terms.length === 0) {
    return []
  }

  return new Promise((resolve, reject) => {
    api.getGroups({search_terms})
      .then(resolve)
      .catch(error =>
        reject(new Error(`Unable to fetch the groups by search terms (${error})`))
      )
  })
}

export async function fetchUsersByNamePrefix(api, prefix) {
  if (prefix.length === 0) {
    return []
  }

  return new Promise((resolve, reject) => {
    api.getUsers({ prefix })
      .then(users => {
        const lowerPrefix = prefix.toLowerCase()
        const withQueryInName = users.filter(
          user => user.name.toLowerCase().indexOf(lowerPrefix) !== -1
        )
        withQueryInName.sort((a, b) => {
          const aValue = a.name.toLowerCase()
          const bValue = b.name.toLowerCase()
          if (aValue.startsWith(lowerPrefix)) {
            return -1
          } else if (bValue.startsWith(lowerPrefix)) {
            return 1
          } else {
            return 0
          }
        })
        resolve(withQueryInName.slice(0, 5))
      })
      .catch(reject)
  })
}
