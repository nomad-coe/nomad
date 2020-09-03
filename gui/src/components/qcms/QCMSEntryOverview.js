import React from 'react'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography } from '@material-ui/core'

export default function QCMSEntryOverview(props) {
  if (!props.data) {
    return <Typography color="error">No metadata available</Typography>
  }
  return (
    <Quantity column>
      <Quantity row>
        <Quantity quantity="formula" label="formula" noWrap {...props} />
        <Quantity quantity="qcms.chemical" label="chemical name" noWrap {...props} />
      </Quantity>
    </Quantity>
  )
}

QCMSEntryOverview.propTypes = {
  data: PropTypes.object.isRequired,
  loading: PropTypes.bool
}
