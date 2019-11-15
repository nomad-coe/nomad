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
  }

  state = {
    preparingDownload: false
  }

  async onDownloadClicked() {
    const {api, query, user, fileName, raiseError} = this.props
    const url = new URL(`${apiBase}/raw/query`)
    url.searchParams.append('strip', 'true')
    Object.keys(query).forEach(key => {url.searchParams.append(key, query[key])})

    if (user) {
      try {
        const token = await api.getSignatureToken()
        url.searchParams.append('signature_token', token)
      } catch(e) {
        this.setState({preparingDownload: false})
        raiseError(e)
      }
    }
    FileSaver.saveAs(url.href, fileName || 'nomad-download.zip')
    this.setState({preparingDownload: false})
  }

  render() {
    const {tooltip, disabled, buttonProps} = this.props
    const {preparingDownload} = this.state

    const props = {
      ...buttonProps,
      disabled: disabled || preparingDownload,
      onClick: () => this.onDownloadClicked()
    }

    return <IconButton {...props}>
      <Tooltip title={tooltip || 'Download'}>
        <DownloadIcon />
      </Tooltip>
    </IconButton>
  }
}

export default compose(withApi(false), withErrors)(DownloadButton)
