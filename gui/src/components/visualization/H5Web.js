
//  TODO: Create a local H5Web component here
import React from 'react'
import PropTypes from 'prop-types'
import { appBase } from '../../config'
import { App, H5GroveProvider } from '@h5web/app'

const H5Web = ({filepath}) => {
  return (
    <H5GroveProvider url={appBase + '/h5grove/'} filepath={filepath} axiosParams={{file: filepath}}>
      <App initialPath={'/entry/optional_parent/optional_child'} />
    </H5GroveProvider>
  )
}
H5Web.propTypes = {
  filepath: PropTypes.string.isRequired
}

export default H5Web
