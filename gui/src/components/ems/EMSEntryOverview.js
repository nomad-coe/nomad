import React from 'react'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Typography, Link } from '@material-ui/core'
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
    const { ems } = data

    if (!ems) {
      return <Typography color="error">No metadata available</Typography>
    }

    const preview_url = ems && ems.preview_url

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
        {data.ems.experiment_summary && <Quantity quantity="ems.experiment_summary" label="summary" {...this.props} />}
        {this.state.previewBroken
          ? data.ems.entry_repository_url && <Quantity label="preview" {...this.props}>
            <Typography noWrap>
              <Link target="external" href={ems.entry_repository_url}>visit this entry on the external database</Link>
            </Typography>
          </Quantity>
          : <Quantity label="preview" {...this.props}>
            <img alt="preview" style={{maxWidth: '100%', height: 'auto'}} src={relative_preview_url} onError={this.handleBrokenPreview}></img>
          </Quantity>}
        <Quantity row>
          <Quantity column>
            <Quantity row>
              <Quantity quantity="formula" label="sample formula" noWrap {...this.props} />
              {data.ems.chemical !== 'unavailable'
                ? <Quantity quantity="ems.chemical" label="sample chemical" noWrap {...this.props} />
                : ''}
            </Quantity>
            <Quantity quantity="ems.method" label="experimental method" noWrap {...this.props} />
            {data.ems.experiment_location && <Quantity quantity="ems.experiment_location" label="experiment location" noWrap {...this.props} />}
            <Quantity label="experiment or experiment publish date" {...this.props}>
              <Typography noWrap>{
                (ems && ems.origin_time && new Date(ems.origin_time).toLocaleDateString()) || 'unavailable'
              }</Typography>
            </Quantity>
            <Quantity label="data source" {...this.props}>
              <Typography noWrap>
                <Link target="external" href={ems.entry_repository_url}>{ems.repository_url}</Link>
              </Typography>
            </Quantity>
          </Quantity>
        </Quantity>
      </Quantity>
    )
  }
}
