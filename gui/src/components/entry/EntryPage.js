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
import React, { useRef } from 'react'
import { Tab, Tabs } from '@material-ui/core'
import { trimEnd } from 'lodash'
import OverviewView from './OverviewView'
import ArchiveEntryView from './ArchiveEntryView'
import ArchiveLogView from './ArchiveLogView'
import BrowseEntryFilesView from './BrowseEntryFilesView'
import { useRouteMatch, useHistory, matchPath, Redirect } from 'react-router-dom'
import { CacheRoute, CacheSwitch } from 'react-router-cache-route'
import EntryContext from './EntryContext'

export const help = `
The *overview* tab gives you an insightful overview about the most prominent
contents found in an entry. You can find more details in the *archive* tab.

The *raw data* tab will show you all files that belong to the entry and offers a download
on individual, or all files. The files can be selected and downloaded. You can also
view the contents of some files directly here on this page.

The *archive* tab shows you the parsed data as a tree
data structure. This view is connected to NOMAD's [metainfo](/metainfo), which acts a schema for
all parsed data.

The *log* tab will show you a log of the entry's processing.
`

const EntryPage = React.memo(() => {
  const history = useHistory()
  const currentPath = history.location.pathname
  const {path, url} = useRouteMatch()
  const urlNoSlash = trimEnd(url, '/')
  const match = matchPath(currentPath, { path: `${path}/:tab?` })
  const {params: {tab = 'overview'}} = match
  const entryId = match?.params?.entryId

  // We use a useRef object to keep track of the current urls of each tab. Switching
  // tabs would go to the previous tab url. This way, the views behind a tab can add
  // state to the URL (e.g. path to section on the ArchiveEntryView).
  const urls = useRef({
    'overview': `${urlNoSlash}/overview`,
    'raw': `${urlNoSlash}/raw`,
    'browse': `${urlNoSlash}/browse`,
    'archive': `${urlNoSlash}/archive`,
    'logs': `${urlNoSlash}/logs`
  })

  const handleChange = (_, value) => {
    urls.current[tab] = currentPath
    history.push(urls.current[value])
  }

  return <EntryContext entryId={entryId}>
    <Tabs
      value={tab}
      onChange={handleChange}
      indicatorColor="primary"
      textColor="primary"
      variant="fullWidth"
    >
      <Tab label="Overview" value="overview" />
      <Tab label="Raw data" value="raw" />
      <Tab label="Processed data" value="archive"/>
      <Tab label="Logs" value="logs"/>
    </Tabs>
    <CacheSwitch>
      <CacheRoute path={`${path}`} exact render={() => <OverviewView/>} />
      <CacheRoute when="always" path={`${path}/raw`} render={() => <BrowseEntryFilesView />} />
      <CacheRoute when="back" path={`${path}/archive`} render={() => <ArchiveEntryView />} />
      <CacheRoute path={`${path}/logs`} render={() => <ArchiveLogView />} />
      <Redirect strict from={`${path}/overview`} to={`${path}`} />
    </CacheSwitch>
  </EntryContext>
})

export default EntryPage
