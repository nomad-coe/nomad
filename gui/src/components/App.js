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
import { Helmet } from "react-helmet"
import { Router, Route } from 'react-router-dom'
import { QueryParamProvider } from 'use-query-params'
import { RecoilRoot } from 'recoil'
import CssBaseline from '@material-ui/core/CssBaseline'
import history from '../history'
import DateFnsUtils from '@date-io/date-fns'
import { MuiPickersUtilsProvider } from '@material-ui/pickers'
import { ui, nomadTheme, keycloakBase, keycloakRealm, keycloakClientId } from '../config'
import Keycloak from 'keycloak-js'
import { ReactKeycloakProvider } from '@react-keycloak/web'
import { MuiThemeProvider } from '@material-ui/core/styles'
import { ErrorSnacks, ErrorBoundary } from './errors'
import Navigation from './nav/Navigation'
import GUIMenu from './GUIMenu'
import { APIProvider, GlobalLoginRequired, onKeycloakEvent } from './api'
import DataStore from './DataStore'
import { UnitProvider } from './units/UnitContext'
import { GlobalMetainfo } from './archive/metainfo'

const keycloak = new Keycloak({
  url: keycloakBase,
  realm: keycloakRealm,
  clientId: keycloakClientId
})

export default function App() {
  return <React.Fragment>
    <Helmet>
      <title>{ui?.theme?.title || 'NOMAD'}</title>
    </Helmet>
    <ReactKeycloakProvider
      authClient={keycloak}
      onEvent={onKeycloakEvent(keycloak)}
      initOptions={{ onLoad: 'check-sso', 'checkLoginIframe': false, promiseType: 'native' }}
      LoadingComponent={<div />}
    >
      <RecoilRoot>
        <APIProvider>
          <MuiPickersUtilsProvider utils={DateFnsUtils}>
            <ErrorSnacks>
              <ErrorBoundary>
                <UnitProvider
                  initialUnitSystems={ui?.unit_systems?.options}
                  initialSelected={ui?.unit_systems?.selected}
                  >
                  <DataStore>
                    <GlobalMetainfo>
                      <Router history={history}>
                        <QueryParamProvider ReactRouterRoute={Route}>
                          <MuiThemeProvider theme={nomadTheme}>
                            <CssBaseline />
                            <GlobalLoginRequired>
                              <Navigation />
                              <GUIMenu/>
                            </GlobalLoginRequired>
                          </MuiThemeProvider>
                        </QueryParamProvider>
                      </Router>
                    </GlobalMetainfo>
                  </DataStore>
                </UnitProvider>
              </ErrorBoundary>
            </ErrorSnacks>
          </MuiPickersUtilsProvider>
        </APIProvider>
      </RecoilRoot>
    </ReactKeycloakProvider>
  </React.Fragment>
}
