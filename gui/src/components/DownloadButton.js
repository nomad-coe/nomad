import React from 'react'
import PropTypes from 'prop-types'
import FileSaver from 'file-saver'
import { withApi } from './api'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { apiBase } from '../config'
import { Tooltip, IconButton, Menu, MenuItem } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'

class DownloadButton extends React.Component {
  static propTypes = {
    /**
     * The query that defines what to download.
     */
    query: PropTypes.object.isRequired,
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
    preparingDownload: false,
    anchorEl: null
  }

  handleClick(event) {
    event.stopPropagation()
    this.setState({ anchorEl: event.currentTarget })
  }

  async handleSelect(choice) {
    const {api, query, user, raiseError} = this.props

    const params = {}
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

    // TODO using split/append is not consistent for all parameters, using append here
    const urlSearchParams = new URLSearchParams()
    Object.keys(params).forEach(key => {
      const value = params[key]
      if (Array.isArray(value)) {
        value.forEach(item => urlSearchParams.append(key, item))
      } else {
        urlSearchParams.append(key, value)
      }
    })
    const url = `${apiBase}/${choice}/${choice === 'archive' ? 'download' : 'query'}?${urlSearchParams.toString()}`
    FileSaver.saveAs(url, `nomad-${choice}-download.zip`)
    this.setState({preparingDownload: false, anchorEl: null})
  }

  handleClose() {
    this.setState({anchorEl: null})
  }

  render() {
    const {tooltip, disabled, buttonProps, dark} = this.props
    const {preparingDownload, anchorEl} = this.state

    const props = {
      ...buttonProps,
      disabled: disabled || preparingDownload,
      onClick: this.handleClick.bind(this)
    }

    return <React.Fragment>
      <IconButton {...props} style={dark ? {color: 'white'} : null}>
        <Tooltip title={tooltip || 'Download'}>
          <DownloadIcon />
        </Tooltip>
      </IconButton>
      <Menu
        anchorEl={anchorEl}
        open={Boolean(anchorEl)}
        onClose={this.handleClose.bind(this)}
      >
        <MenuItem onClick={() => this.handleSelect('raw')}>Raw uploaded files</MenuItem>
        <MenuItem onClick={() => this.handleSelect('archive')}>NOMAD Archive files</MenuItem>
      </Menu>
    </React.Fragment>
  }
}

export default compose(withApi(false), withErrors)(DownloadButton)
