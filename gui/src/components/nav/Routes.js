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
import { Route, useLocation } from 'react-router-dom'
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
// import SearchPage, {help as searchHelp} from '../search/SearchPage'
import NewSearchPage, {help as newSearchHelp} from '../search/NewSearchPage'
import UploadPage, {help as uploadHelp} from '../uploads/UploadPage'
import UserdataPage, {help as userdataHelp} from '../UserdataPage'
import { ErrorBoundary } from '../errors'

function createEntryRoute(props) {
  return ({
    path: '/entry',
    title: 'Entry',
    help: {
      title: 'The entry page',
      content: entryHelp
    },
    ...props,
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
  })
}

const routeSpecs = [
  {
    path: '/faq',
    exact: true,
    title: 'Frequently Asked Questions',
    component: FAQ
  },
  {
    path: '/search',
    exact: true,
    title: 'Search and Download Data',
    help: {
      title: 'How to find and download data',
      content: newSearchHelp
    },
    navPath: 'explore/search',
    component: NewSearchPage
  },
  {
    path: '/userdata',
    exact: true,
    title: 'Manage Your Data',
    help: {
      title: 'How to manage your data',
      content: userdataHelp
    },
    navPath: 'publish/userdata',
    component: UserdataPage,
    routes: [
      createEntryRoute({
        navPath: 'publish/userdata',
        breadCrumbs: [
          {
            title: 'Your Data',
            path: '/userdata'
          }
        ]
      })
    ]
  },
  createEntryRoute({
    navPath: 'explore/search',
    breadCrumbs: [
      {
        title: 'Search',
        path: '/search'
      }
    ]
  }),
  {
    path: '/dataset',
    title: 'Dataset',
    navPath: 'explore/search',
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
  {
    path: '/uploads',
    exact: true,
    title: 'Upload and Publish Data',
    help: {
      title: 'How to upload data',
      content: uploadHelp
    },
    navPath: 'publish/uploads',
    component: UploadPage,
    routes: [
      createEntryRoute({
        navPath: 'publish/uploads',
        breadCrumbs: [
          {
            title: 'Uploads',
            path: '/uploads'
          }
        ]
      })
    ]
  },
  {
    path: '/metainfo',
    title: 'The NOMAD Meta Info',
    help: {
      title: 'About the NOMAD meta-info',
      content: metainfoHelp
    },
    navPath: 'analyze/metainfo',
    component: MetainfoPage
  },
  {
    path: '/aitoolkit',
    title: 'Artificial Intelligence Toolkit',
    navPath: 'analyze/aitoolkit',
    component: AIToolkitPage
  },
  {
    exact: true,
    path: '/apis',
    title: 'APIs',
    navPath: 'analyze/apis',
    component: APIs
  },
  {
    exact: true,
    path: '/',
    title: 'About, Documentation, Getting Help',
    navPath: 'about/info',
    component: About
  }
]

function flattenRouteSpecs(routeSpecs, parent, results) {
  results = results || []
  parent = parent || {}
  routeSpecs.forEach(route => {
    const flatRoute = {
      ...parent,
      component: null,
      exact: false,
      ...route,
      path: parent.path ? `${parent.path}/${route.path.replace(/^\/+/, '')}` : route.path,
      routes: undefined
    }

    if (flatRoute.component) {
      results.push(flatRoute)
    }

    if (route.routes) {
      flattenRouteSpecs(route.routes, flatRoute, results)
    }
  })
  return results
}

export const routes = flattenRouteSpecs(routeSpecs)
routes.sort((a, b) => (a.path > b.path) ? -1 : 1)

export default function Routes() {
  return <React.Fragment>
    {routes.map(route => {
      const {path, exact} = route
      const children = childProps => childProps.match && <route.component {...childProps} />
      return <ErrorBoundary key={path}>
        <Route exact={exact} path={path}
          // eslint-disable-next-line react/no-children-prop
          children={children}
        />
      </ErrorBoundary>
    })}
  </React.Fragment>
}

export function useRoute() {
  const {pathname, search} = useLocation()
  routes.forEach(route => {
    (route.breadCrumbs || []).forEach(breadCrumb => {
      if (breadCrumb.path.startsWith(pathname)) {
        breadCrumb.path = pathname + (search || '')
      }
    })
  })
  const route = routes.find(route => pathname.startsWith(route.path))
  return route
}
