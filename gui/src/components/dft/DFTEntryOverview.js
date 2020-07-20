import React from 'react'
import PropTypes from 'prop-types'
import { Typography, Button, makeStyles, Tooltip } from '@material-ui/core'
import Quantity from '../Quantity'
import _ from 'lodash'
import {appBase} from '../../config'

const useStyles = makeStyles(theme => ({
  actions: {
    marginTop: theme.spacing(1),
    textAlign: 'right',
    margin: -theme.spacing(1)
  }
}))

export default function DFTEntryOverview(props) {
  const classes = useStyles()
  const {data} = props
  if (!data.dft) {
    return <Typography color="error">No metadata available</Typography>
  }

  const material_name = entry => entry.encyclopedia.material.material_name

  return <div>
    <Quantity column>
      <Quantity row>
        <Quantity quantity="formula" label='formula' noWrap {...props} />
        <Quantity quantity={material_name} label='material name' noWrap {...props} />
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
    {data.encyclopedia && data.encyclopedia.material &&
      <div className={classes.actions}>
        <Tooltip title="Show the material of this entry in the NOMAD Encyclopedia.">
          <Button color="primary" href={`${appBase}/encyclopedia/#/material/${data.encyclopedia.material.material_id}`}>
            material
          </Button>
        </Tooltip>
      </div>
    }
  </div>
}

DFTEntryOverview.propTypes = {
  data: PropTypes.object.isRequired
}
