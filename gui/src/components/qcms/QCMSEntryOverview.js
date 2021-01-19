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
import { Typography } from '@material-ui/core'
import { apiBase } from '../../config'

/**
 * Shows an informative overview about the selected entry.
 */
export default function QCMSEntryOverview({repo, uploadId, calcId}) {
  const [previewBroken, setPreviewBroken] = useState(false)
  const { qcms } = repo
  if (!qcms) {
    return <Typography color="error">No metadata available</Typography>
  }
  const dirname = repo.mainfile.substring(0, repo.mainfile.lastIndexOf('/'))
  const relative_preview_url = `${apiBase}/raw/${repo.upload_id}/${dirname}/circuit.png`

  return (
    <DefaultEntryOverview repo={repo} uploadId={uploadId} calcId={calcId}>
      <Quantity column>
        <Quantity row>
          <Quantity quantity="formula" label="formula" noWrap data={repo}/>
          <Quantity quantity="qcms.chemical" label="chemical name" noWrap data={repo}/>
        </Quantity>
        <Quantity row>
          <Quantity quantity="qcms.quantum_computer_system" label="system" noWrap data={repo}/>
          <Quantity quantity="qcms.quantum_computing_libraries" label="libraries" noWrap data={repo}/>
        </Quantity>
        <Quantity label="compute time">
          <Typography noWrap>{new Date(qcms.computation_datetime).toLocaleDateString()}</Typography>
        </Quantity>
        {!previewBroken &&
         <Quantity label="circuit" data={repo}>
           <img alt="circuit" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={() => setPreviewBroken(true)}></img>
         </Quantity>}
      </Quantity>
    </DefaultEntryOverview>
  )
}

QCMSEntryOverview.propTypes = {
  repo: PropTypes.object.isRequired,
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
