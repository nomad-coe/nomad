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
import { Typography } from '@material-ui/core'
import { apiBase } from '../../config'

export default function QCMSEntryDetails({data}) {
  const [previewBroken, setPreviewBroken] = useState(false)

  if (!data) {
    return <Typography color="error">No metadata available</Typography>
  }

  let dirname = data.mainfile.substring(0, data.mainfile.lastIndexOf('/'))
  if (dirname !== '') dirname += '/'
  const relative_preview_url = `${apiBase}/raw/${data.upload_id}/${dirname}circuit.png`

  return (
    <Quantity column>
      <Quantity row>
        <Quantity quantity="formula" label="formula" noWrap data={data}/>
        <Quantity quantity="qcms.chemical" label="chemical name" noWrap data={data}/>
      </Quantity>
      <Quantity row>
        <Quantity quantity="qcms.quantum_computer_system" label="system" noWrap data={data}/>
        <Quantity quantity="qcms.quantum_computing_libraries" label="libraries" noWrap data={data}/>
      </Quantity>
      <Quantity label="compute time">
        <Typography noWrap>{new Date(data.qcms.computation_datetime).toLocaleDateString()}</Typography>
      </Quantity>
      {!previewBroken &&
        <Quantity label="circuit" data={data}>
          <img alt="circuit" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={() => setPreviewBroken(true)}></img>
        </Quantity>}
    </Quantity>
  )
}

QCMSEntryDetails.propTypes = {
  data: PropTypes.object.isRequired
}
