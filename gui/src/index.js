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
import 'react-app-polyfill/ie11'
import 'react-app-polyfill/stable'
import React from 'react'
import ReactDOM from 'react-dom'
import './index.css'
import App from './components/App'
import { Router, Route } from 'react-router-dom'
import { QueryParamProvider } from 'use-query-params'
import history from './history'
import PiwikReactRouter from 'piwik-react-router'
import { matomoEnabled, matomoUrl, matomoSiteId, keycloakBase, keycloakRealm, keycloakClientId } from './config'
import Keycloak from 'keycloak-js'
import { KeycloakProvider } from 'react-keycloak'
import * as serviceWorker from './serviceWorker'

export const matomo = matomoEnabled ? PiwikReactRouter({
  url: matomoUrl,
  siteId: matomoSiteId,
  clientTrackerName: 'stat.js',
  serverTrackerName: 'stat'
}) : []

const keycloak = Keycloak({
  url: keycloakBase,
  realm: keycloakRealm,
  clientId: keycloakClientId
})

// matomo.push('requireConsent')

ReactDOM.render(
  <KeycloakProvider keycloak={keycloak} initConfig={{onLoad: 'check-sso', 'checkLoginIframe': false}} LoadingComponent={<div />}>
    <Router history={matomoEnabled ? matomo.connectToHistory(history) : history}>
      <QueryParamProvider ReactRouterRoute={Route}>
        <App />
      </QueryParamProvider>
    </Router>
  </KeycloakProvider>, document.getElementById('root'))

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
serviceWorker.unregister()
