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
import React, { } from 'react'
import PropTypes from 'prop-types'
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'
import { Trajectories } from './Trajectory'

/**
 * Shows a summary of thermodynamic properties for the given entry data.
 */
const ThermodynamicProperties = React.memo(({index, archive}) => {
  return <PropertyGrid>
    <PropertyItem xs={12} title="Trajectory" height="auto">
      <Trajectories index={index} archive={archive}/>
    </PropertyItem>
  </PropertyGrid>
})

ThermodynamicProperties.propTypes = {
  index: PropTypes.any,
  archive: PropTypes.any
}

export default ThermodynamicProperties
