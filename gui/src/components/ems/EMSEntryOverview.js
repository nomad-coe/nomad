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

  state = {
    previewBroken: false
  }

  constructor(props) {
    super(props)
    this.handleBrokenPreview = this.handleBrokenPreview.bind(this)
  }

  handleBrokenPreview(event) {
    this.setState({previewBroken: true})
  }

  render() {
    const { data } = this.props
    const { ems: { preview_url } } = data

    let relative_preview_url = null
    if (!preview_url) {
      relative_preview_url = 'broken'
    } else if (preview_url.indexOf('http://') === 0 || preview_url.indexOf('https://') === 0) {
      relative_preview_url = preview_url
    } else {
      const dirname = data.mainfile.substring(0, data.mainfile.lastIndexOf('/'))
      relative_preview_url = `${apiBase}/raw/${data.upload_id}/${dirname}/${preview_url}`
    }

    return (
      <Quantity column>
        <Quantity quantity="ems.experiment_summary" label="summary" {...this.props} />
        <Quantity row>
          <Quantity column>
            <Quantity row>
              <Quantity quantity="formula" label="sample formula" noWrap {...this.props} />
              {data.ems.chemical !== 'unavailable'
                ? <Quantity quantity="ems.chemical" label="sample chemical" noWrap {...this.props} />
                : ''}
            </Quantity>
            <Quantity quantity="ems.method" label="experimental method" noWrap {...this.props} />
            <Quantity quantity="ems.experiment_location" label="experiment location" noWrap {...this.props} />
            <Quantity label="experiment time" {...this.props}>
              <Typography noWrap>{
                data.ems.experiment_time !== 'unavailable' ? new Date(data.ems.experiment_time * 1000).toLocaleString() : 'unavailable'
              }</Typography>
            </Quantity>
            <Quantity label="data" {...this.props}>
              <Typography noWrap>
                <a target="external" href={data.ems.repository_url}>{data.ems.repository_name}</a>
              </Typography>
            </Quantity>
          </Quantity>
          {this.state.previewBroken
            ? data.ems.entry_repository_url && <Quantity label="preview" {...this.props}>
              <Typography noWrap>
                <a target="external" href={data.ems.entry_repository_url}>visit this entry on the external database</a>
              </Typography>
            </Quantity>
            : <Quantity label="preview" {...this.props}>
              <img alt="preview" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={this.handleBrokenPreview}></img>
            </Quantity>}
        </Quantity>
      </Quantity>
    )
  }
}
