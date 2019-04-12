import React from 'react'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography } from '@material-ui/core'

export default class EMSEntryOverview extends React.Component {
  static propTypes = {
    data: PropTypes.object.isRequired,
    loading: PropTypes.bool
  }

  render() {
    const { data } = this.props

    return (
      <Quantity row>
        <Quantity column>
          <Quantity row>
            <Quantity quantity="formula" label="sample formula" noWrap {...this.props} />
            <Quantity quantity="chemical" label="sample chemical" noWrap {...this.props} />
          </Quantity>
          <Quantity quantity="method" label="experimental method" noWrap {...this.props} />
          <Quantity quantity="experiment_location" label="experiment location" noWrap {...this.props} />
          <Quantity label="experiment time" {...this.props}>
            <Typography noWrap>
              {new Date(data.experiment_time * 1000).toLocaleString()}
            </Typography>
          </Quantity>
          <Quantity label="data" {...this.props}>
            <Typography noWrap>
              <a href={data.repository_url}>{data.repository_name}</a>
            </Typography>
          </Quantity>
        </Quantity>
        <Quantity label="preview" {...this.props}>
          <img alt="preview" style={{maxWidth: '100%', height: 'auto'}} src={data.preview_url}></img>
        </Quantity>
      </Quantity>
    )
  }
}
