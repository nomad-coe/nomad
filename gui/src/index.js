import 'react-app-polyfill/ie11'
import 'react-app-polyfill/stable'
import React from 'react'
import ReactDOM from 'react-dom'
import './index.css'
import App from './components/App'
import registerServiceWorker from './registerServiceWorker'
import { Router, Route } from 'react-router-dom'
import { QueryParamProvider } from 'use-query-params'
import history from './history'
import PiwikReactRouter from 'piwik-react-router'
import { sendTrackingData, matomoUrl, matomoSiteId, keycloakBase, keycloakRealm, keycloakClientId } from './config'
import Keycloak from 'keycloak-js'
import { KeycloakProvider } from 'react-keycloak'

const matomo = sendTrackingData ? PiwikReactRouter({
  url: matomoUrl,
  siteId: matomoSiteId
}) : null

const keycloak = Keycloak({
  url: keycloakBase,
  realm: keycloakRealm,
  clientId: keycloakClientId
})

ReactDOM.render(
  <KeycloakProvider keycloak={keycloak} initConfig={{onLoad: 'check-sso'}} LoadingComponent={<div />}>
    <Router history={sendTrackingData ? matomo.connectToHistory(history) : history}>
      <QueryParamProvider ReactRouterRoute={Route}>
        <App />
      </QueryParamProvider>
    </Router>
  </KeycloakProvider>, document.getElementById('root'))
registerServiceWorker()
