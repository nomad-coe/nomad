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
import React, { useState } from 'react'
import DefaultEntryOverview from '../entry/DefaultEntryOverview'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography, Link } from '@material-ui/core'
import { apiBase } from '../../config'

/**
 * Shows an informative overview about the selected entry.
 */
export default function EMSEntryOverview({repo, uploadId, calcId}) {
  const [previewBroken, setPreviewBroken] = useState(false)
  const handleBrokenPreview = (event) => {
    setPreviewBroken(true)
  }
  const { ems } = repo
  if (!ems) {
    return <Typography color="error">No metadata available</Typography>
  }

  const preview_url = ems && ems.preview_url
  let relative_preview_url = null
  if (!preview_url) {
    relative_preview_url = 'broken'
  } else if (preview_url.indexOf('http://') === 0 || preview_url.indexOf('https://') === 0) {
    relative_preview_url = preview_url
  } else {
    const dirname = repo.mainfile.substring(0, repo.mainfile.lastIndexOf('/'))
    relative_preview_url = `${apiBase}/raw/${repo.upload_id}/${dirname}/${preview_url}`
  }

  return (
    <DefaultEntryOverview repo={repo} uploadId={uploadId} calcId={calcId}>
      <Quantity column>
        {ems.experiment_summary && <Quantity quantity="ems.experiment_summary" label="summary" data={repo}/>}
        {previewBroken
          ? ems.entry_repository_url && <Quantity label="preview" data={repo}>
            <Typography noWrap>
              <Link target="external" href={ems.entry_repository_url}>visit this entry on the external database</Link>
            </Typography>
          </Quantity>
          : <Quantity label="preview" data={repo}>
            <img alt="preview" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={handleBrokenPreview}></img>
          </Quantity>}
        <Quantity row>
          <Quantity column>
            <Quantity row>
              <Quantity quantity="formula" label="sample formula" noWrap data={repo}/>
              {ems.chemical !== 'unavailable'
                ? <Quantity quantity="ems.chemical" label="sample chemical" noWrap data={repo}/>
                : ''}
            </Quantity>
            <Quantity quantity="ems.method" label="experimental method" noWrap data={repo}/>
            {ems.experiment_location && <Quantity quantity="ems.experiment_location" label="experiment location" noWrap data={repo}/>}
            <Quantity label="experiment or experiment publish date" data={repo}>
              <Typography noWrap>{
                (ems && ems.origin_time && new Date(ems.origin_time).toLocaleDateString()) || 'unavailable'
              }</Typography>
            </Quantity>
            <Quantity label="data source" data={repo}>
              <Typography noWrap>
                <Link target="external" href={ems.entry_repository_url}>{ems.repository_url}</Link>
              </Typography>
            </Quantity>
          </Quantity>
        </Quantity>
      </Quantity>
    </DefaultEntryOverview>
  )
}

EMSEntryOverview.propTypes = {
  repo: PropTypes.object.isRequired,
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
