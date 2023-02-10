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
import { useRouteMatch, useHistory, matchPath, Redirect } from 'react-router-dom'
import { CacheRoute, CacheSwitch } from 'react-router-cache-route'
import UploadOverview from './UploadOverview'
import UploadFilesView from './UploadFilesView'
import UploadPageContext from './UploadPageContext'
import UploadProcessingStatus from './UploadProcessingStatus'
import PropTypes from 'prop-types'

const UploadPage = React.memo((props) => {
  const history = useHistory()
  const { uploadId: propsUploadId } = props
  const currentPath = history.location.pathname
  const {path, url} = useRouteMatch()
  const urlNoSlash = trimEnd(url, '/')
  const match = matchPath(currentPath, { path: `${path}/:tab?` })
  const {params: {tab = 'overview'}} = match
  const uploadId = propsUploadId || match?.params?.uploadId

  // We use a useRef object to keep track of the current urls of each tab. Switching
  // tabs would go to the previous tab url. This way, the views behind a tab can add
  // state to the URL (e.g. path to section on the ArchiveEntryView).
  const urls = useRef({
    'overview': `${urlNoSlash}/overview`,
    'files': `${urlNoSlash}/files`
  })

  const handleChange = (_, value) => {
    urls.current[tab] = currentPath
    history.push(urls.current[value])
  }
  // TODO I removed the UploadOverview view from the cache, because
  // processing in the UploadFile (Browser) view do not update the
  // upload context and do not refresh the UploadOverview data,
  // Test
  return <UploadPageContext uploadId={uploadId}>
    <UploadProcessingStatus />
    <Tabs
      value={tab}
      onChange={handleChange}
      indicatorColor="primary"
      textColor="primary"
      variant="fullWidth"
    >
      <Tab label="Overview" value="overview" />
      <Tab label="Files" value="files" />
    </Tabs>
    {tab === 'overview' && (
      <UploadOverview/>
    )}
    <CacheSwitch>
      <CacheRoute when="always" path={`${path}/files`} render={() => <UploadFilesView />} />
      <Redirect strict from={`${path}/overview`} to={`${path}`} />
    </CacheSwitch>
  </UploadPageContext>
})
UploadPage.propTypes = {
  uploadId: PropTypes.string
}

export default UploadPage
