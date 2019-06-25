import React from 'react'
import PropTypes from 'prop-types'
import { Typography } from '@material-ui/core'
import Quantity from '../Quantity'

export default class DFTEntryOverview extends React.Component {
  static propTypes = {
    data: PropTypes.object.isRequired,
    loading: PropTypes.bool
  }

  render() {
    const { data } = this.props

    return (
      <Quantity column>
        <Quantity row>
          <Quantity quantity="code_name" label='dft code' noWrap {...this.props} />
          <Quantity quantity="code_version" label='dft code version' noWrap {...this.props} />
        </Quantity>
        <Quantity row>
          <Quantity quantity="basis_set" label='basis set' noWrap {...this.props} />
          <Quantity quantity="xc_functional" label='xc functional' noWrap {...this.props} />
        </Quantity>
        <Quantity row>
          <Quantity quantity="system" label='system type' noWrap {...this.props} />
          <Quantity quantity="crystal_system" label='crystal system' noWrap {...this.props} />
          <Quantity quantity='spacegroup_symbol' label="spacegroup" noWrap {...this.props}>
            <Typography noWrap>
              {data.spacegroup_symbol} ({data.spacegroup})
            </Typography>
          </Quantity>
        </Quantity>
      </Quantity>
    )
  }
}
