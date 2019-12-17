import React from 'react'
import PropTypes from 'prop-types'
import FileSaver from 'file-saver'
import { withApi } from './api'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { apiBase } from '../config'
import { Tooltip, IconButton } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'

class DownloadButton extends React.Component {
  static propTypes = {
    /**
     * The query that defines what to download.
     */
    query: PropTypes.object.isRequired,
    /**
     * A suggestion for the download filename.
     */
    fileName: PropTypes.string,
    /**
     * A tooltip for the button
     */
    tooltip: PropTypes.string,
    /**
     * Whether the button is disabled
     */
    disabled: PropTypes.bool,
    /**
     * Properties forwarded to the button.
     */
    buttonProps: PropTypes.object,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    raiseError: PropTypes.func.isRequired,
    dark: PropTypes.bool
  }

  state = {
    preparingDownload: false
  }

  async onDownloadClicked(event) {
    event.stopPropagation()

    const {api, query, user, fileName, raiseError} = this.props

    const params = {
      strip: true
    }
    Object.keys(query).forEach(key => { params[key] = query[key] })

    if (user) {
      try {
        const token = await api.getSignatureToken()
        params['signature_token'] = token
      } catch (e) {
        this.setState({preparingDownload: false})
        raiseError(e)
      }
    }
    FileSaver.saveAs(`${apiBase}/raw/query?${new URLSearchParams(params).toString()}`, fileName || 'nomad-download.zip')
    this.setState({preparingDownload: false})
  }

  render() {
    const {tooltip, disabled, buttonProps, dark} = this.props
    const {preparingDownload} = this.state

    const props = {
      ...buttonProps,
      disabled: disabled || preparingDownload,
      onClick: this.onDownloadClicked.bind(this)
    }

    return <IconButton {...props} style={dark ? {color: 'white'} : null}>
      <Tooltip title={tooltip || 'Download'}>
        <DownloadIcon />
      </Tooltip>
    </IconButton>
  }
}

export default compose(withApi(false), withErrors)(DownloadButton)
