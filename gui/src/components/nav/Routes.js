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
import { Route } from 'react-router-dom'
import About from '../About'
import APIs from '../APIs'
import AIToolkitPage from '../aitoolkit/AIToolkitPage'
import { MetainfoPage, help as metainfoHelp } from '../archive/MetainfoBrowser'
import ResolveDOI from '../dataset/ResolveDOI'
import DatasetPage from '../DatasetPage'
import EntryPage, {help as entryHelp} from '../entry/EntryPage'
import EntryQuery from '../entry/EntryQuery'
import ResolvePID from '../entry/ResolvePID'
import FAQ from '../FAQ'
import SearchPage, {help as searchHelp} from '../search/SearchPage'
import UploadPage, {help as uploadHelp} from '../uploads/UploadPage'
import UserdataPage, {help as userdataHelp} from '../UserdataPage'

export const routes = {
  'faq': {
    path: '/faq',
    exact: true,
    appBarTitle: 'Frequently Asked Questions',
    component: FAQ
  },
  'search': {
    path: '/search',
    exact: true,
    appBarTitle: 'Find and Download Data',
    appBarHelp: {
      title: 'How to find and download data',
      content: searchHelp
    },
    navPath: 'explore/search',
    component: SearchPage
  },
  'userdata': {
    path: '/userdata',
    exact: true,
    appBarTitle: 'Manage Your Data',
    appBarHelp: {
      title: 'How to manage your data',
      content: userdataHelp
    },
    navPath: 'publish/userdata',
    component: UserdataPage
  },
  'entry': {
    path: '/entry',
    appBarTitle: 'Entry',
    appBarHelp: {
      title: 'The entry page',
      content: entryHelp
    },
    defaultNavPath: 'explore/search',
    routes: [
      {
        path: '/id',
        component: EntryPage
      },
      {
        path: '/query',
        exact: true,
        component: EntryQuery
      },
      {
        path: '/pid',
        component: ResolvePID
      }
    ]
  },
  'dataset': {
    path: '/dataset',
    appBarTitle: 'Dataset',
    defaultNavPath: 'explore/search',
    routes: [
      {
        path: '/id',
        component: DatasetPage
      },
      {
        path: '/doi',
        component: ResolveDOI
      }
    ]
  },
  'uploads': {
    path: '/uploads',
    exact: true,
    appBarTitle: 'Upload and Publish Data',
    appBarHelp: {
      title: 'How to upload data',
      content: uploadHelp
    },
    navPath: 'publish/uploads',
    component: UploadPage
  },
  'metainfo': {
    path: '/metainfo',
    appBarTitle: 'The NOMAD Meta Info',
    appBarHelp: {
      title: 'About the NOMAD meta-info',
      content: metainfoHelp
    },
    navPath: 'analyze/metainfo',
    component: MetainfoPage
  },
  'aitoolkit': {
    path: '/aitoolkit',
    appBarTitle: 'Artificial Intelligence Toolkit',
    navPath: 'analyze/aitoolkit',
    component: AIToolkitPage
  },
  'apis': {
    exact: true,
    path: '/apis',
    appBarTitle: 'APIs',
    navPath: 'analyze/apis',
    component: APIs
  },
  'about': {
    exact: true,
    path: '/',
    appBarTitle: 'About, Documentation, Getting Help',
    navPath: 'about/info',
    component: About
  }
}

export const allRoutes = Object.keys(routes).map(key => {
  const route = routes[key]
  return route.routes
    ? route.routes.map(subRoute => ({
      ...route,
      ...subRoute,
      path: `${route.path}/${subRoute.path.replace(/^\/+/, '')}`,
      routes: undefined
    }))
    : route
}).flat()

export default function Routes() {
  return <React.Fragment>
    {Object.keys(allRoutes).map(routeKey => {
      const route = allRoutes[routeKey]
      const { path, exact } = route
      const children = childProps => childProps.match && <route.component {...childProps} />
      return <Route key={routeKey} exact={exact} path={path}
        // eslint-disable-next-line react/no-children-prop
        children={children}
      />
    })}
  </React.Fragment>
}
