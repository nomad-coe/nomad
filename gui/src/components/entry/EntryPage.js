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
import PropTypes from 'prop-types'
import { Tab, Tabs, Box } from '@material-ui/core'
import OverviewView from './OverviewView'
import ArchiveEntryView from './ArchiveEntryView'
import ArchiveLogView from './ArchiveLogView'
import RawFileView from './RawFileView'
import { useRouteMatch, useHistory, matchPath, Route } from 'react-router-dom'

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

const TabRoutes = React.memo(function Routes({match}) {
  const {params: {uploadId, entryId, tab = 'overview'}} = match
  const props = {entryId: entryId, uploadId: uploadId}

  if (!entryId) {
    return ''
  }

  return <React.Fragment>
    <Box display={tab === 'overview' ? 'block' : 'none'}><OverviewView {...props}/></Box>
    <Box display={tab === 'raw' ? 'block' : 'none'}><RawFileView {...props}/></Box>
    <Box display={tab === 'archive' ? 'block' : 'none'}><ArchiveEntryView {...props}/></Box>
    <Box display={tab === 'logs' ? 'block' : 'none'}><ArchiveLogView {...props}/></Box>
  </React.Fragment>
})

TabRoutes.propTypes = {
  match: PropTypes.object.isRequired
}

const EntryPage = React.memo(function EntryPage() {
  const history = useHistory()
  const currentPath = history.location.pathname
  const {path, url} = useRouteMatch()

  const match = matchPath(currentPath, { path: `${path}/:tab?` })
  const {params: {tab = 'overview'}} = match

  // We use a useRef object to keep track of the current urls of each tab. Switching
  // tabs would go to the previous tab url. This way, the views behind a tab can add
  // state to the URL (e.g. path to section on the ArchiveEntryView).
  const urls = useRef({
    'overview': `${url}/overview`,
    'raw': `${url}/raw`,
    'archive': `${url}/archive`,
    'logs': `${url}/logs`
  })

  const handleChange = (_, value) => {
    urls.current[tab] = currentPath
    history.push(urls.current[value])
  }

  return <React.Fragment>
    <Tabs
      value={tab} onChange={handleChange}
      indicatorColor="primary" textColor="primary" variant="fullWidth"
    >
      <Tab label="Overview" value="overview" />
      <Tab label="Raw data" value="raw" />
      <Tab label="Archive" value="archive"/>
      <Tab label="Logs" value="logs"/>
    </Tabs>
    <Route path={`${path}/:tab?`} render={(props) => <TabRoutes {...props}/>} />
  </React.Fragment>
})

export default EntryPage
