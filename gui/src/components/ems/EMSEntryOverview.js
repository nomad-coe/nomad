import React from 'react'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography } from '@material-ui/core'
import { apiBase } from '../../config'

export default class EMSEntryOverview extends React.Component {
  static propTypes = {
    data: PropTypes.object.isRequired,
    loading: PropTypes.bool
  }

  render() {
    const { data } = this.props
    const { preview_url } = data

    let relative_preview_url = null
    if (!preview_url) {
      relative_preview_url = 'broken'
    } else if (preview_url.indexOf('http://') === 0 || preview_url.indexOf('https://') === 0) {
      relative_preview_url = preview_url
    } else {
      const dirname = data.mainfile.substring(0, data.mainfile.lastIndexOf('/') + 1)
      relative_preview_url = `${apiBase}/raw/${data.upload_id}/${dirname ? dirname + '/' : ''}${preview_url}`
    }

    return (
      <Quantity column>
        <Quantity quantity="experiment_summary" label="summary" {...this.props} />
        <Quantity row>
          <Quantity column>
            <Quantity row>
              <Quantity quantity="formula" label="sample formula" noWrap {...this.props} />
              {data.chemical !== 'unavailable'
                ? <Quantity quantity="chemical" label="sample chemical" noWrap {...this.props} />
                : ''}
            </Quantity>
            <Quantity quantity="method" label="experimental method" noWrap {...this.props} />
            <Quantity quantity="experiment_location" label="experiment location" noWrap {...this.props} />
            <Quantity label="experiment time" {...this.props}>
              <Typography noWrap>{
                data.experiment_time !== 'unavailable' ? new Date(data.experiment_time * 1000).toLocaleString() : 'unavailable'
              }</Typography>
            </Quantity>
            <Quantity label="data" {...this.props}>
              <Typography noWrap>
                <a href={data.repository_url}>{data.repository_name}</a>
              </Typography>
            </Quantity>
          </Quantity>
          <Quantity label="preview" {...this.props}>
            <img alt="preview" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url}></img>
          </Quantity>
        </Quantity>
      </Quantity>
    )
  }
}
