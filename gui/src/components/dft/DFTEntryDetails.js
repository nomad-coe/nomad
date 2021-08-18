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
import { Typography, Tooltip, Link } from '@material-ui/core'
import Quantity from '../Quantity'
import _ from 'lodash'
import {appBase, encyclopediaEnabled, normalizeDisplayValue} from '../../config'

export default function DFTEntryDetails({data}) {
  if (!data?.dft) {
    return <Typography color="error">No metadata available</Typography>
  }

  const material_name = entry => {
    let name
    try {
      name = entry.encyclopedia.material.material_name
    } catch {}
    name = name || 'unnamed'

    if (encyclopediaEnabled && data?.encyclopedia?.material?.material_id && data.published && !data.with_embargo) {
      const url = `${appBase}/encyclopedia/#/material/${data.encyclopedia.material.material_id}`
      return (
        <Tooltip title="Show the material of this entry in the NOMAD Encyclopedia.">
          <Link href={url}>{name}</Link>
        </Tooltip>
      )
    } else {
      return name
    }
  }

  return <div>
    <Quantity column>
      <Quantity row>
        <Quantity quantity="formula" label='formula' noWrap data={data}/>
        <Quantity quantity={material_name} label='material' noWrap data={data}/>
      </Quantity>
      <Quantity row>
        <Quantity quantity="dft.code_name" label='dft code' noWrap data={data}/>
        <Quantity quantity="dft.code_version" label='dft code version' noWrap data={data}/>
      </Quantity>
      <Quantity row>
        <Quantity quantity="dft.basis_set" label='basis set' noWrap data={data}/>
        <Quantity quantity="dft.xc_functional" label='xc functional' noWrap data={data}/>
      </Quantity>
      <Quantity row>
        <Quantity quantity="dft.system" label='system type' noWrap data={data}/>
        <Quantity quantity="dft.crystal_system" label='crystal system' noWrap data={data}/>
        <Quantity quantity="dft.spacegroup_symbol" label="spacegroup" noWrap data={data}>
          <Typography noWrap>
            {normalizeDisplayValue(_.get(data, 'dft.spacegroup_symbol'))} ({normalizeDisplayValue(_.get(data, 'dft.spacegroup'))})
          </Typography>
        </Quantity>
      </Quantity>
    </Quantity>
  </div>
}

DFTEntryDetails.propTypes = {
  data: PropTypes.object.isRequired
}
