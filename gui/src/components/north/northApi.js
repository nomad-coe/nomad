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

import axios from 'axios'
import { northBase } from '../../config'

export class NorthApi {
  constructor(api, base) {
    this.api = api
    this.apiKey = null
    this.axios = axios.create({
      baseURL: `${northBase}/hub/api${base ? '/' + base : ''}`
    })
  }

  async getApiKey() {
    if (this.apiKey) {
      return this.apiKey
    }

    if (!this.api.keycloak.token) {
      throw Error('User does not have a token. This should not happen')
    }

    const config = {
      headers: {
        accept: 'application/json',
        Authorization: `Bearer ${this.api.keycloak.token}`
      }
    }

    const response = await this.axios.post('tokens', null, config)
    this.apiKey = response.data.token

    return this.apiKey
  }

  async get(path, query, config) {
    const method = (path, body, config) => this.axios.get(path, config)
    return this.doHttpRequest(method, path, null, {params: query, ...config})
  }

  async post(path, body, config) {
    const method = (path, body, config) => this.axios.post(path, body, config)
    return this.doHttpRequest(method, path, body, config)
  }

  async put(path, body, config) {
    const method = (path, body, config) => this.axios.put(path, body, config)
    return this.doHttpRequest(method, path, body, config)
  }

  async delete(path, config) {
    const method = (path, body, config) => this.axios.delete(path, config)
    return this.doHttpRequest(method, path, null, config)
  }

  async doHttpRequest(method, path, body, config) {
    const apiKey = await this.getApiKey()
    config = config || {}
    config.params = config.params || {}
    config.headers = config.headers || {
      accept: 'application/json',
      Authorization: `token ${apiKey}`
    }
    try {
      return method(path, body, config)
    } catch (errors) {
      if (config.noHandleErrors) {
        throw errors
      }
    }
  }
}
