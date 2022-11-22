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
import FileBrowser from '../archive/FileBrowser'
import Page from '../Page'
import { useUploadPageContext } from './UploadPageContext'
import { createUploadUrl } from '../../utils'
import { Typography } from '@material-ui/core'

const UploadFilesView = React.memo(function UploadFilesView() {
  const {deploymentUrl, uploadId, error, hasUpload} = useUploadPageContext()

  if (!hasUpload) {
    return <Page limitedWidth>
      {(error ? <Typography>{error.apiMessage || error.message || 'Failed to load'}</Typography> : <Typography>loading ...</Typography>)}
    </Page>
  }

  return <Page>
    <FileBrowser
      uploadUrl={createUploadUrl(deploymentUrl, uploadId, '')}
      rootTitle="Upload files"
    />
  </Page>
})
UploadFilesView.propTypes = {
  uploadId: PropTypes.string
}
export default UploadFilesView
