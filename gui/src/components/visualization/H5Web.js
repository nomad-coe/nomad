
//  TODO: Create a local H5Web component here
import React from 'react'
import { appBase } from '../../config';
import { App, H5GroveProvider } from '@h5web/app';

const H5Web = (props) => {
  return(
    <H5GroveProvider url={appBase+"/h5grove/"} filepath={props.filepath} axiosParams={{file: props.filepath}}>
        <App initialPath={"/entry/optional_parent/optional_child"} />
    </H5GroveProvider>
  )
}

export default H5Web