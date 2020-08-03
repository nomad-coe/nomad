import React from 'react'
import PropTypes from 'prop-types'
import { Typography, Tooltip, Link } from '@material-ui/core'
import Quantity from '../Quantity'
import _ from 'lodash'
import {appBase, encyclopediaEnabled} from '../../config'

export default function DFTEntryOverview(props) {
  const {data} = props
  if (!data.dft) {
    return <Typography color="error">No metadata available</Typography>
  }

  const material_name = entry => {
    let name
    try {
      name = entry.encyclopedia.material.material_name
    } catch {}
    name = name || 'unnamed'

    if (encyclopediaEnabled && data.encyclopedia && data.encyclopedia.material && data.published && !data.with_embargo) {
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
        <Quantity quantity="formula" label='formula' noWrap {...props} />
        <Quantity quantity={material_name} label='material' noWrap {...props} />
      </Quantity>
      <Quantity row>
        <Quantity quantity="dft.code_name" label='dft code' noWrap {...props} />
        <Quantity quantity="dft.code_version" label='dft code version' noWrap {...props} />
      </Quantity>
      <Quantity row>
        <Quantity quantity="dft.basis_set" label='basis set' noWrap {...props} />
        <Quantity quantity="dft.xc_functional" label='xc functional' noWrap {...props} />
      </Quantity>
      <Quantity row>
        <Quantity quantity="dft.system" label='system type' noWrap {...props} />
        <Quantity quantity="dft.crystal_system" label='crystal system' noWrap {...props} />
        <Quantity quantity="dft.spacegroup_symbol" label="spacegroup" noWrap {...props}>
          <Typography noWrap>
            {_.get(data, 'dft.spacegroup_symbol')} ({_.get(data, 'dft.spacegroup')})
          </Typography>
        </Quantity>
      </Quantity>
    </Quantity>
  </div>
}

DFTEntryOverview.propTypes = {
  data: PropTypes.object.isRequired
}
