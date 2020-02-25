import React from 'react'
import PropTypes from 'prop-types'
import { Typography } from '@material-ui/core'
import Quantity from '../Quantity'
import _ from 'lodash'

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
          <Quantity quantity="formula" label='formula' noWrap {...this.props} />
        </Quantity>
        <Quantity row>
          <Quantity quantity="dft.code_name" label='dft code' noWrap {...this.props} />
          <Quantity quantity="dft.code_version" label='dft code version' noWrap {...this.props} />
        </Quantity>
        <Quantity row>
          <Quantity quantity="dft.basis_set" label='basis set' noWrap {...this.props} />
          <Quantity quantity="dft.xc_functional" label='xc functional' noWrap {...this.props} />
        </Quantity>
        <Quantity row>
          <Quantity quantity="dft.system" label='system type' noWrap {...this.props} />
          <Quantity quantity="dft.crystal_system" label='crystal system' noWrap {...this.props} />
          <Quantity quantity="dft.spacegroup_symbol" label="spacegroup" noWrap {...this.props}>
            <Typography noWrap>
              {_.get(data, 'dft.spacegroup_symbol')} ({_.get(data, 'dft.spacegroup')})
            </Typography>
          </Quantity>
        </Quantity>
      </Quantity>
    )
  }
}
