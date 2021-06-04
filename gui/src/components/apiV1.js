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
import React from 'react'
import PropTypes from 'prop-types'
import { apiBase } from '../config'
import { makeStyles, Typography } from '@material-ui/core'
import LoginLogout from './LoginLogout'
import { useKeycloak } from 'react-keycloak'
import axios from 'axios'

export class DoesNotExist extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'DoesNotExist'
  }
}

export class NotAuthorized extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'NotAuthorized'
  }
}

export class ApiError extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'CannotReachApi'
  }
}

function handleApiError(e) {
  if (e.name === 'CannotReachApi' || e.name === 'NotAuthorized' || e.name === 'DoesNotExist') {
    throw e
  }

  let error = null
  if (e.response) {
    const body = e.response.body
    const message = (body && (body.message || body.description)) || e.response.statusText
    const errorMessage = `${message} (${e.response.status})`
    if (e.response.status === 404) {
      error = new DoesNotExist(errorMessage)
    } else if (e.response.status === 401) {
      error = new NotAuthorized(errorMessage)
    } else if (e.response.status === 502) {
      error = new ApiError(errorMessage)
    } else {
      error = new Error(errorMessage)
    }
    error.status = e.response.status
    error.apiMessage = message
  } else {
    if (e.message === 'Failed to fetch') {
      error = new ApiError(e.message)
      error.status = 400
    } else {
      const errorMessage = e.status ? `${e} (${e.status})` : '' + e
      error = new Error(errorMessage)
    }
  }
  throw error
}

class Api {
  constructor(keycloak) {
    this.keycloak = keycloak
    this.axios = axios.create({
      baseURL: `${apiBase}/v1`
    })

    this.loadingHandler = []
    this.loading = 0

    this.onFinishLoading = () => {
      this.loading--
      this.loadingHandler.forEach(handler => handler(this.loading))
    }
    this.onStartLoading = () => {
      this.loading++
      this.loadingHandler.forEach(handler => handler(this.loading))
    }
  }

  /**
   * Fetches the up-to-date authorization headers. Performs a token update if
   * the current token has expired.
   * @returns Object containing the authorization header.
   */
  async authHeaders() {
    if (this.keycloak.token) {
      return new Promise((resolve, reject) => {
        this.keycloak.updateToken()
          .success(() => {
            resolve({
              headers: {
                'Authorization': `Bearer ${this.keycloak.token}`
              }
            })
          })
          .error(() => {
            reject(new ApiError())
          })
      })
    } else {
      return {}
    }
  }

  /**
   * Returns the entry information that is stored in the search index.
   *
   * @param {string} entryId
   * @returns Object containing the search index contents.
   */
  async entry(entryId) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const entry = await this.axios.get(`/entries/${entryId}`, auth)
      return entry.data
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }

  /**
   * Returns section_results from the archive corresponding to the given entry.
   * Some large quantities which are not required by the GUI are filtered out.
   * All references within the section are resolved by the server before
   * sending.
   *
   * @param {string} entryId
   * @returns Object containing section_results
   */
  async results(entryId) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const entry = await this.axios.post(
        `/entries/${entryId}/archive/query`,
        {
          required: {
            'resolve-inplace': false,
            results: {
              material: '*',
              method: '*',
              properties: {
                structures: '*',
                electronic: 'include-resolved',
                vibrational: 'include-resolved',
                // We do not want to include the resolved trajectory in the
                // response, so we explicitly specify which parts we want
                geometry_optimization: {
                  energies: 'include-resolved'
                }
              }
            }
          }
        },
        auth
      )
      return parse(entry).data.data.archive
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }

  /**
   * Returns a list of suggestions for the given metainfo quantity.
   *
   * @param {string} quantity
   * @returns List of suggested values as dictionaries containing both the
   * suggestion and its weight. The items are ordered by weight.
   */
  async suggestions(quantities, input) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const suggestions = await this.axios.post(
        `/suggestions`,
        {
          input: input,
          quantities: quantities
        },
        auth
      )
      return suggestions.data
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }

  /**
   * Return the raw file metadata for a given entry.
   * @param {string} entryId
   * @returns Object containing the raw file metadata.
   */
  async getRawFileListFromCalc(entryId) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const entry = await this.axios.get(`/entries/${entryId}/raw`, auth)
      return entry.data
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }

  /**
   * Executes the given entry query
   * @param {object} query contains the query
   * @returns Object containing the raw file metadata.
   */
  async queryEntry(search) {
    this.onStartLoading()
    const auth = await this.authHeaders()
    try {
      const result = await this.axios.post(
        '/entries/query',
        {
          exclude: ['atoms', 'only_atoms', 'files', 'dft.quantities', 'dft.optimade', 'dft.labels', 'dft.geometries'],
          ...search
        },
        auth
      )
      return result.data
    } catch (errors) {
      handleApiError(errors)
    } finally {
      this.onFinishLoading()
    }
  }
}

/**
 * Custom JSON parse function that can handle NaNs that can be created by the
 * python JSON serializer.
 */
function parse(result) {
  if (typeof result === 'string') {
    try {
      return JSON.parse(result)
    } catch (e) {
      try {
        return JSON.parse(result.replace(/\bNaN\b/g, '"NaN"'))
      } catch (e) {
        return result
      }
    }
  } else {
    return result
  }
}

let api = null

export function useApi() {
  const [keycloak] = useKeycloak()

  if (!api || api.keycloak !== keycloak) {
    api = new Api(keycloak)
  }

  return api
}

const useLoginRequiredStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(2),
    display: 'flex',
    alignItems: 'center',
    '& p': {
      marginRight: theme.spacing(1)
    }
  }
}))

export function LoginRequired({message, children}) {
  const classes = useLoginRequiredStyles()
  const api = useApi()
  if (api.keycloak.authenticated) {
    return <React.Fragment>
      {children}
    </React.Fragment>
  } else {
    return <div className={classes.root}>
      <Typography>
        {message || 'You have to login to use this functionality.'}
      </Typography>
      <LoginLogout color="primary" />
    </div>
  }
}
LoginRequired.propTypes = {
  message: PropTypes.string,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
}
