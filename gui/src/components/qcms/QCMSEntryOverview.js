import React, { useState } from 'react'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography } from '@material-ui/core'
import { apiBase } from '../../config'

export default function QCMSEntryOverview(props) {
  const [previewBroken, setPreviewBroken] = useState(false)

  if (!props.data) {
    return <Typography color="error">No metadata available</Typography>
  }

  const {data} = props
  const dirname = data.mainfile.substring(0, data.mainfile.lastIndexOf('/'))
  const relative_preview_url = `${apiBase}/raw/${data.upload_id}/${dirname}/circuit.png`

  return (
    <Quantity column>
      <Quantity row>
        <Quantity quantity="formula" label="formula" noWrap {...props} />
        <Quantity quantity="qcms.chemical" label="chemical name" noWrap {...props} />
      </Quantity>
      <Quantity row>
        <Quantity quantity="qcms.quantum_computer_system" label="system" noWrap {...props} />
        <Quantity quantity="qcms.quantum_computing_libraries" label="libraries" noWrap {...props} />
      </Quantity>
      <Quantity label="compute time">
        <Typography noWrap>{new Date(data.qcms.computation_datetime).toLocaleDateString()}</Typography>
      </Quantity>
      {!previewBroken &&
        <Quantity label="circuit" {...props}>
          <img alt="circuit" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={() => setPreviewBroken(true)}></img>
        </Quantity>}
    </Quantity>
  )
}

QCMSEntryOverview.propTypes = {
  data: PropTypes.object.isRequired,
  loading: PropTypes.bool
}
