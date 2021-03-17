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
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography, Link } from '@material-ui/core'
import { apiBase } from '../../config'

export default function EMSEntryDetails({data}) {
  const [previewBroken, setPreviewBroken] = useState(false)
  const handleBrokenPreview = (event) => {
    setPreviewBroken(true)
  }
  const { ems } = data
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
    const dirname = data.mainfile.substring(0, data.mainfile.lastIndexOf('/'))
    relative_preview_url = `${apiBase}/raw/${data.upload_id}/${dirname}/${preview_url}`
  }

  return (
    <Quantity column>
      {ems.experiment_summary && <Quantity quantity="ems.experiment_summary" label="summary" data={data}/>}
      {previewBroken
        ? ems.entry_repository_url && <Quantity label="preview" data={data}>
          <Typography noWrap>
            <Link target="external" href={ems.entry_repository_url}>visit this entry on the external database</Link>
          </Typography>
        </Quantity>
        : <Quantity label="preview" data={data}>
          <img alt="preview" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={handleBrokenPreview}></img>
        </Quantity>}
      <Quantity row>
        <Quantity column>
          <Quantity row>
            <Quantity quantity="formula" label="sample formula" noWrap data={data}/>
            {ems.chemical !== 'unavailable'
              ? <Quantity quantity="ems.chemical" label="sample chemical" noWrap data={data}/>
              : ''}
          </Quantity>
          <Quantity quantity="ems.method" label="experimental method" noWrap data={data}/>
          {ems.experiment_location && <Quantity quantity="ems.experiment_location" label="experiment location" noWrap data={data}/>}
          <Quantity label="experiment or experiment publish date" data={data}>
            <Typography noWrap>{
              (ems && ems.origin_time && new Date(ems.origin_time).toLocaleDateString()) || 'unavailable'
            }</Typography>
          </Quantity>
          <Quantity label="data source" data={data}>
            <Typography noWrap>
              <Link target="external" href={ems.entry_repository_url}>{ems.repository_url}</Link>
            </Typography>
          </Quantity>
        </Quantity>
      </Quantity>
    </Quantity>
  )
}

EMSEntryDetails.propTypes = {
  data: PropTypes.object.isRequired
}
