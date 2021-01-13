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
import { Tab, Tabs, Box } from '@material-ui/core'
import OverviewView from './OverviewView'
import ArchiveEntryView from './ArchiveEntryView'
import ArchiveLogView from './ArchiveLogView'
import RepoEntryView from './RepoEntryView'
import KeepState from '../KeepState'
import { guiBase } from '../../config'
import { useRouteMatch, useHistory, Route } from 'react-router-dom'

export const help = `
The *raw files* tab, will show you all files that belong to the entry and offers a download
on individual, or all files. The files can be selected and downloaded. You can also
view the contents of some files directly here on this page.

The *archive* tab, shows you the parsed data as a tree
data structure. This view is connected to NOMAD's [meta-info](${guiBase}/metainfo), which acts a schema for
all parsed data.

The *log* tab, will show you a log of the entry's processing.
`

export function EntryPageContent({children, fixed}) {
  const props = fixed ? {maxWidth: 1024} : {}
  return <Box padding={3} margin="auto" {...props}>
    {children}
  </Box>
}
EntryPageContent.propTypes = ({
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  fixed: PropTypes.bool
})

export default function EntryPage() {
  const history = useHistory()
  const { path, url } = useRouteMatch()

  return (
    <Route
      path={`${path}/:uploadId?/:calcId?/:tab?`}
      render={({match: {params: {uploadId, calcId, tab = 'raw'}}}) => {
        if (calcId && uploadId) {
          const calcProps = { calcId: calcId, uploadId: uploadId }
          return (
            <React.Fragment>
              <Tabs
                value={tab || 'overview'}
                onChange={(_, value) => history.push(`${url}/${uploadId}/${calcId}/${value}`)}
                indicatorColor="primary"
                textColor="primary"
                variant="fullWidth"
              >
                <Tab label="Overview" value="overview" />
                <Tab label="Raw data" value="raw" />
                <Tab label="Archive" value="archive"/>
                <Tab label="Logs" value="logs"/>
              </Tabs>

              <KeepState visible={tab === 'overview' || tab === undefined} render={props => <OverviewView {...props} />} {...calcProps} />
              <KeepState visible={tab === 'raw'} render={props => <RepoEntryView {...props} />} {...calcProps} />
              <KeepState visible={tab === 'archive'} render={props => <ArchiveEntryView {...props} />} {...calcProps} />
              <KeepState visible={tab === 'logs'} render={props => <ArchiveLogView {...props} />} {...calcProps} />
            </React.Fragment>
          )
        } else {
          return ''
        }
      }}
    />
  )
}
